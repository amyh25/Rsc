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


#' run_monocle_trajectory
#' 
#' Wrapper for the monocle pipeline: 
#' so %>% as.cell_data_set() %>% cluster_cells() %>% learn_graph() %>% order_cells()
#' 
#' @param so
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