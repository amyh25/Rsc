#' make_pseudobulk
#' 
#' Makes pseudobulk from Seurat object
#' 
#' @param so Seurat object
#' @param split_var as string, name of variable splitting pseudobulks (e.g. cluster)
#' @param sample_var as string, name of variable identifying samples
#' @param vars list of strings, name of variables of relevant conditions
#' @return list of pseudobulks (i.e. matrices of aggregated counts)
#' 
#' 
#' @export

make_pseudobulk <- function(so, split_var = NULL, sample_var, vars) {
  
  sce <- seurat_to_sce(so, split_var, sample_var, vars)
  
  if (is.null(split_var)) {
    
    groups <- SummarizedExperiment::colData(sce)[, c(sample_var, vars)]
    pb <- Matrix.utils::aggregate.Matrix(
      Matrix::t(SingleCellExperiment::counts(sce)), 
      groupings = groups, fun = "sum"
    )
    return(pb)
        
  } else {

    groups <- map(sce, ~SummarizedExperiment::colData(..1)[, vars])
    pb_list <- map2(
      sce, groups, 
      ~Matrix.utils::aggregate.Matrix(
        Matrix::t(SingleCellExperiment::counts(..1)), 
        groupings = ..2, fun = "sum"
      )
    ) %>% set_names(names(sce))
    
    return(pb_list)
    
  }
}

#' qc_sce
#' 
#' QC for sce
#' 
#' @param sce_obj SingleCellExperiment object
#' @param sample_var as string, name of variable identifying samples
#' @param threshold min counts per gene threshold (default: 10)
#' @return sce object
#' @export

qc_sce <- function(sce_obj, sample_var, threshold = 10) {
  sample_names <- unique(sce_obj[[sample_var]]) %>% set_names()
  n_cells <- as.numeric(table(sce_obj[[sample_var]]))
  m <- base::match(sample_names, sce_obj[[sample_var]])
  sample_level_metadata <- 
    SummarizedExperiment::colData(sce_obj)[m, ] %>% 
    data.frame(., n_cells, row.names = NULL)
  
  sce_cell <- scater::perCellQCMetrics(sce_obj)
  sce_obj$is_outlier <- scater::isOutlier(
    metric = sce_cell$total, 
    nmads = 2, type = "both", log = TRUE
  )
  sce_obj <- sce_obj[, !sce_obj$is_outlier]
  sce_obj <- sce_obj[Matrix::rowSums(SingleCellExperiment::counts(sce_obj) > 1) >= threshold, ]
  
  return(sce_obj)
}

#' seurat_to_sce
#' 
#' Turns seurat objects into sce objects
#' 
#' @param so Seurat object
#' @param split_var as string, name of variable splitting pseudobulks (e.g. cluster)
#' @param sample_var as string, name of variable identifying samples
#' @param vars list of strings, name of variables of relevant condition
#' @return sce object
#' @export

seurat_to_sce <- function(so, split_var=NULL, sample_var, vars) {
  
  seurat_counts <- so@assays$RNA@counts
  seurat_metadata <- so@meta.data %>% 
    dplyr::select(sample_var, vars)
  
  if (is.null(split_var)) {
    
    sce_obj <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = seurat_counts), 
      colData = seurat_metadata
    )
    sce_obj <- qc_sce(sce_obj, sample_var)
    
    return(sce_obj)
    
  } else {
    
    seurat_metadata <- so@meta.data %>% 
      dplyr::select(vars, sample_var, split_var)
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = seurat_counts), 
      colData = seurat_metadata
    )
    
    sce_list <- unique(seurat_metadata[[split_var]]) %>% 
      purrr::map(~base::subset(sce, , get(split_var) == ..1)) %>% 
      set_names(unique(seurat_metadata[[split_var]]))
    
    # qc sce cluster list
    sce_list <- purrr::map(sce_list, qc_sce, sample_var)
    
    return(sce_list)
    
  }
  
}
