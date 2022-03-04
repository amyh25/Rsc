#' scp_export_clustering_file
#' 
#' Create clustering file for export to the Single Cell Portal
#' https://singlecell.broadinstitute.org/single_cell
#' 
#' @param so Seurat object
#' @param str_subset_regex_group regex string to select columns from Seurat object metadata (default: "RNA_snn_res.")
#' @param str_subset_regex_numeric regex string to select columns from Seurat object metadata (default: "_RNA")
#' @param reduction_str string specifying name of reduction. Should match the parameter reduction.name from Seurat::RunUMAP (default: "umap")
#' @param output_dir path to output directory
#' @param filename export filename (default: "clustering_scp")
#' @return Seurat object with identity set to given ident
#' @export

scp_export_clustering_file <- function(so, 
                                       str_subset_regex_group = "RNA_snn_res.", 
                                       str_subset_regex_numeric = "_RNA", 
                                       reduction_str = "umap", 
                                       output_dir, 
                                       filename = "clustering_scp") {
  
  umap_coords <- so@reductions[[reduction_str]]@cell.embeddings %>% 
    as_tibble(rownames = "cell")
  metadata_df <- so@meta.data %>% 
    as_tibble(rownames = "cell") %>% 
    left_join(umap_coords, by = "cell")
  
  other_cols_to_select_group <- str_subset(colnames(metadata_df), str_subset_regex_group)
  other_cols_to_select_num <- str_subset(colnames(metadata_df), str_subset_regex_numeric)
  clustering_cols <- c(
    "NAME", "X", "Y", 
    other_cols_to_select_group, 
    other_cols_to_select_num
  )
  
  clustering_types_vec <- c(
    "TYPE", "numeric", "numeric", 
    rep("group", (length(other_cols_to_select_group))), 
    rep("numeric", (length(other_cols_to_select_num)))
  ) %>% set_names(clustering_cols)
  
  clustering_types_line <- paste0(clustering_types_vec, collapse = "\t")
  
  clustering_scp <- metadata_df %>% 
    mutate(NAME = cell, X = UMAP_1, Y = UMAP_2) %>% 
    select(names(clustering_types_vec))
  write_tsv(clustering_scp, file.path(output_dir, "tmp.txt"))
  
  f_clustering_raw <- file(file.path(output_dir, "tmp.txt"), open = "r")
  lines <- readLines(f_clustering_raw)
  close(f_clustering_raw)
  
  new_lines <- c(lines[1], clustering_types_line, lines[2:length(lines)])
  f_clustering <- file(file.path(output_dir, paste0(filename, ".txt")), open = "w")
  writeLines(new_lines, f_clustering)
  close(f_clustering)
  
  file.remove(file.path(output_dir, "tmp.txt"))
  
  print("done!")
  
}

