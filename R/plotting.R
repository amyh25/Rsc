#' make_galaxy_plot
#' 
#' Makes galaxy plot
#' 
#' @param .data data frame with umap coordinates
#' @param sample_n number of cells to downsample to (optional; default: 800)
#' @param umap1 string denoting UMAP_1 (optional; default: UMAP_1)
#' @param umap2 string denoting UMAP_2 (optional; default: UMAP_2)
#' @return 
#' @export

make_galaxy_plot <- function(.data, 
                             sample_n = 800, 
                             umap1 = "UMAP_1", umap2 = "UMAP_2") {
  data_sampled <- slice_sample(data, n = sample_n)
  ggplot(data) +
    aes_string(umap1, umap2) + 
    stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE) +
    scale_fill_viridis(option = "magma") +
    coord_cartesian(expand = FALSE, xlim = c(min(data[[umap1]]), max(data[[umap1]])),
                    ylim = c(min(data[[umap2]]), max(data[[umap2]]))) + 
    geom_point(shape = '.', col = 'white', data = data_sampled)
}

#' plot_binary_on_umap
#' 
#' Plots a binary variable on UMAP
#' 
#' @param metadata_df data frame with umap coordinates and binary variable
#' @param binary_var binary variable, must be TRUE/FALSE
#' @param umap1 string denoting UMAP_1 (optional; default: UMAP_1)
#' @param umap2 string denoting UMAP_2 (optional; default: UMAP_2)
#' @return 
#' @export

plot_binary_on_umap <- function(metadata_df, binary_var, 
                                umap1 = "UMAP_1", umap2 = "UMAP_2") {
  binary_var <- enquo(binary_var)
  count_df <- metadata_df %>% 
    group_by(!!binary_var, orig.ident) %>% 
    summarise(ncells = n(), .groups = "drop")
  metadata_df %>% 
    ggplot() + 
    aes_string(umap1, umap2) + 
    aes(color = !!binary_var) + 
    geom_point(size = 0.5, stroke = 0.1, alpha = 0.5) + 
    geom_text(aes(x = -8, y = 10 + 2*as.integer(!!binary_var), label = ncells), 
              data = count_df) + 
    coord_fixed() + 
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey"))
}

#' plot_volcano_skeleton
#' 
#' Plots skeleton of volcano plot
#' 
#' @param .data data frame with p values and fold change
#' @param fc name of fold change variable
#' @param pval name of p-val variable
#' @param remove_inf_fc (optional; default: TRUE) remove zero and infinite fold change values
#' @param pt.size (optional; default: 0.5) size of points
#' @param fc_thresh (optional; default: 1) fold change threshold to draw lines at
#' @return ggplot of volcano plot
#' @export

plot_volcano_skeleton <- function(.data, fc, pval, 
                                  remove_inf_fc = TRUE, 
                                  pt.size = 0.5, 
                                  fc_thresh = 1) {
  fc <- enquo(fc)
  pval <- enquo(pval)
  
  if (remove_inf_fc) {
    .data <- .data %>% 
      filter(!is.infinite(!!fc) & !(!!fc == 0)) 
  }
  
  .data %>% 
    ggplot() + 
    aes(log2(!!fc), -log10(!!pval)) + 
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), color = "grey", size = 0.2) +
    geom_hline(yintercept = -log10(0.05), color = "red", size = 0.2) +
    geom_hline(yintercept = -log10(0.25), color = "pink", size = 0.2) +
    geom_point(size = pt.size, stroke = 0.1) + 
    xlab("Fold Change") + 
    ylab(TeX("-log_{10} (p-val)"))
}
