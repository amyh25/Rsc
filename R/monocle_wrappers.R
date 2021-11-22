#' run_monocle_trajectory
#' 
#' Wrapper for the monocle pipeline: 
#' so %>% as.cell_data_set() %>% cluster_cells() %>% learn_graph() %>% order_cells()
#' 
#' @param so seurat object
#' @return cell data set object
#' @export

run_monocle_trajectory <- function(
  so, reduction_method = "UMAP", 
  use_partition = FALSE, close_loop = FALSE
) {
  so %>% 
    as.cell_data_set() %>% 
    cluster_cells(reduction_method = reduction_method) %>% 
    learn_graph(use_parition = use_partition, 
                close_loop = close_loop) %>% 
    order_cells(reduction_method = reduction_method)
}

#' add_pseudotime_metadata
#' 
#' Adds pseudotime as metadata to seurat object
#' 
#' @param so seurat object
#' @param cds cell data set object
#' @return seurat object with pseudotime metadata
#' @export 

add_pseudotime_metadata <- function(so, cds) {
  as.data.frame(pseudotime(cds))
}