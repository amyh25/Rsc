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
