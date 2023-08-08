#' set_idents
#' 
#' Sets identity of Seurat object
#' Equivalent to Ident(so) <- "ident"
#' 
#' @param so Seurat object
#' @param ident Identity given as string
#' @return Seurat object with identity set to given ident
#' @export

set_idents <- function(so, ident) {
  Idents(so) <- ident
  return(so)
}

#' set_assay
#' 
#' Sets default assay of Seurat object
#' Equivalent to DefaultAssay(so) <- "assay"
#' 
#' @param so Seurat object
#' @param assay Assay given as string
#' @return Seurat object with default assay set to given assay
#' @export

set_assay <- function(so, assay) {
  DefaultAssay(so) <- assay
  return(so)
}

#' add_mito_and_ribo
#' 
#' Adds percentage of mitochondrial reads and ribosomal reads for each cell
#' 
#' @param so Seurat object
#' @param prefix percent string prefix (optional; default: "percent.")
#' @return Seurat object with percentage mitochondrial reads and percent ribosomal reads
#' @export

add_mito_and_ribo <- function(seurat, prefix = "percent.") {
  seurat[[paste0(prefix, "mt")]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
  seurat[[paste0(prefix, "ribo")]] <- PercentageFeatureSet(seurat, pattern = "(^Rpl|^Rps)")
  return(seurat)
}

#' get_gene_expr_from_so
#' 
#' Gets tibble of genes expression for select genes from Seurat object
#' 
#' @param so Seurat object
#' @param gene_vec character vector of gene names
#' @param assay assay to pull data from (optional; default:"RNA")
#' @param skip_missing bool, whether or not to skip missing genes. If true, will error
#' @return tibble of gene expression per cell
#' @export

get_gene_expr_from_so <- function(so, gene_vec, assay = "RNA", skip_missing = FALSE) {
  if (sum(gene_vec %in% rownames(so)) == length(gene_vec)) {
    missing_genes <- gene_vec[!(gene_vec %in% rownames(so))]
    if (skip_missing) {
      message(paste0(missing_genes, " not found. Skipping..."))
      gene_vec <- gene_vec[gene_vec %in% rownames(so)]
    } else {
      stop(paste0(missing_genes, " not found. Please remove. "))
    }
  }
  so@assays[[assay]]@data %>% 
    .[rownames(.) %in% gene_vec,] %>% 
    as.matrix() %>% 
    t() %>% 
    as_tibble(rownames = "cell")
}


#' pivot_FindAllMarkers_wider
#' 
#' Takes output of `FindAllMarkers`, filters down to significant genes, 
#' then builds a wide table of top markers for each ident specified
#' 
#' @param markers Seurat object
#' @param pval pval var (default: p_val_adj)
#' @param pval_thresh p value filtering threshold (default: 0.05)
#' @param log2FC log2 fold change var (default: avg_log2FC)
#' @param ident identity var (default: cluster)
#' @return tibble of top markers per ident, all genes meeting p-value threshold, ordered by top log2FC
#' @export

pivot_FindAllMarkers_wider <- function(markers, 
                                       pval = p_val_adj, 
                                       pval_thresh = 0.05, 
                                       log2FC = avg_log2FC, 
                                       ident = cluster) {
  pval <- enquo(pval)
  log2FC <- enquo(avg_log2FC)
  ident <- enquo(ident)
  markers %>% 
    filter(!!pval < pval_thresh) %>% 
    arrange(desc(!!log2FC)) %>% 
    group_by(!!ident) %>% 
    mutate(rank = row_number()) %>% 
    arrange(!!ident) %>% 
    pivot_wider(id_cols = rank, 
                names_from = !!ident, 
                values_from = gene)
}

