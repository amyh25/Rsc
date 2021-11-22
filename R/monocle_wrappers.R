#' run_monocle_trajectory
#' 
#' Wrapper for the monocle pipeline: 
#' so %>% as.cell_data_set() %>% cluster_cells() %>% learn_graph() %>% order_cells()
#' 
#' @param so seurat object
#' @return cell data set object
#' @export

run_monocle_trajectory <- function(
  so, reduction_method = "UMAP"
) {
  so %>% 
    as.cell_data_set %>% 
    cluster_cells(reduction_method = reduction_method) %>% 
    learn_graph() %>% 
    order_cells(reduction_method = reduction_method)
}