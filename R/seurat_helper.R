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
#' @return tibble of gene expression per cell
#' @export

get_gene_expr_from_so <- function(so, gene_vec, assay = "RNA") {
  so@assays[[assay]]@data %>% 
    .[rownames(.) %in% gene_vec,] %>% 
    as.matrix() %>% 
    t() %>% 
    as_tibble(rownames = "cell")
}

