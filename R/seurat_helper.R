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