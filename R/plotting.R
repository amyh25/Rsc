#' pca_gplot_scale
#' 
#' Plots PCA given a matrix with column names of variables
#' 
#' Function given to me by Deb Sen, credit goes to Kevin Bi
#' 
#' @param matrix Data matrix with column names
#' @param sep Sep character in sample name
#' @param nameindex Part of name that you want to plot
#' @param dim1 PC component on x-axis
#' @param dim2 PC component on y-axis
#' @return 
#' @export

pca_gplot_scale <- function(matrix, sep, nameindex, dim1, dim2){
  library(ggplot2)
  pca_forplot <-  prcomp(t(matrix),center = TRUE, scale = T)
  percent_variation <- (pca_forplot$sdev)^2 / sum(pca_forplot$sdev^2) *100
  batch <-  unlist(lapply(strsplit(colnames(matrix),sep),function(x)x[[nameindex]]))
  xvar <- percent_variation[dim1]
  yvar <- percent_variation[dim2]
  q <- qplot(pca_forplot$x[,dim1], pca_forplot$x[,dim2],colour=batch, 
             xlab = paste("PC", dim1, " (", format(xvar, digits = 4), "% var)"), 
             ylab = paste("PC", dim2, " (", format(yvar, digits = 4), "% var)"), 
             size = I(3))
  q + scale_color_discrete(name = "")
}

#' plot_gene_str
#' 
#' Plots genes, with points ordered by expression of gene value
#' 
#' @param gene_str gene name, as string
#' @param genes_df data frame with at least two columns, one with gene values and one with obs IDs to join to umap_coords
#' @param umap_coords data frame with at least three columns, two with UMAP coords and one with obs IDs to join to genes_df
#' @param UMAP_1 string denoting UMAP_1 (optional; default: UMAP_1)
#' @param UMAP_2 string denoting UMAP_2 (optional; default: UMAP_2)
#' @param stroke_size stroke size (optional; default: 0.1)
#' @param pt_size point size (optional; default: 0.5)
#' @return 
#' @export

plot_gene_str <- function(gene_str, genes_df, umap_coords, 
                          join_col_str = "cell", 
                          UMAP_1 = "UMAP_1", UMAP_2 = "UMAP_2", 
                          stroke_size = 0.1, pt_size = 0.5) {
  genes_df %>% 
    left_join(umap_coords, by = join_col_str) %>% 
    arrange(get(gene_str)) %>% 
    ggplot() + 
    aes_string(UMAP_1, UMAP_2) + 
    aes(color = get(gene_str)) + 
    geom_point(stroke = stroke_size, size = pt_size) + 
    coord_fixed() + 
    xlab("UMAP dim. 1") + ylab("UMAP dim. 2") + 
    labs(color = "Expr.") + 
    ggtitle(gene_str) + 
    scale_color_gradient(low = "grey", high = "blue")
  
}

#' plot_gene
#' 
#' Plots genes, with points ordered by expression of gene value. 
#' If you'd like to use the gene name instead, use `plot_gene_str`
#' 
#' @param gene gene name, passed in directly
#' @param genes_df data frame with at least two columns, one with gene values and one with obs IDs to join to umap_coords
#' @param umap_coords data frame with at least three columns, two with UMAP coords and one with obs IDs to join to genes_df
#' @param UMAP_1 string denoting UMAP_1 (optional; default: UMAP_1)
#' @param UMAP_2 string denoting UMAP_2 (optional; default: UMAP_2)
#' @param stroke_size stroke size (optional; default: 0.1)
#' @param pt_size point size (optional; default: 0.5)
#' @return 
#' @export

plot_gene <- function(gene, genes_df, umap_coords, join_col_str = "cell", 
                      UMAP_1 = "UMAP_1", UMAP_2 = "UMAP_2", 
                      stroke_size = 0.1, pt_size = 0.5) {
  gene <- enquo(gene)
  genes_df %>% 
    left_join(umap_coords, by = join_col_str) %>% 
    arrange(!!gene) %>% 
    ggplot() + 
    aes_string(UMAP_1, UMAP_2) + 
    aes(color = !!gene) + 
    geom_point(stroke = stroke_size, size = pt_size) + 
    coord_fixed() + 
    xlab("UMAP dim. 1") + ylab("UMAP dim. 2") + 
    ggtitle(rlang::as_name(gene)) + 
    scale_color_gradient(low = "grey", high = "blue")
}

#' make_umap_skeleton
#' 
#' Creates umap skeleton from data frame
#' 
#' @param .data data frame with umap coordinates
#' @param UMAP_1 string denoting UMAP_1 (optional; default: UMAP_1)
#' @param UMAP_2 string denoting UMAP_2 (optional; default: UMAP_2)
#' @param stroke_size stroke size (optional; default: 0.1)
#' @param pt_size point size (optional; default: 0.5)
#' @return 
#' @export

make_umap_skeleton <- function(metadata_df, 
                               UMAP_1 = "UMAP_1", UMAP_2 = "UMAP_2", 
                               stroke_size = 0.1, pt_size = 0.5) {
  
  metadata_df %>% 
    ggplot() + 
    aes_string(UMAP_1, UMAP_2) + 
    geom_point(stroke = stroke_size, size = pt_size) + 
    xlab("UMAP dim. 1") + ylab("UMAP dim. 2") + 
    coord_fixed() + 
    guides(color = guide_legend(override.aes = list(size = 3)))
  
}


#' make_galaxy_plot
#' 
#' Makes galaxy plot
#' 
#' @param metadata_df data frame with umap coordinates
#' @param sample_n number of cells to downsample to (optional; default: 800)
#' @param umap1 string denoting UMAP_1 (optional; default: UMAP_1)
#' @param umap2 string denoting UMAP_2 (optional; default: UMAP_2)
#' @return 
#' @export

make_galaxy_plot <- function(metadata_df, 
                             sample_n = 800, 
                             umap1 = "UMAP_1", umap2 = "UMAP_2") {
  if (sample_n != 0)
    data_sampled <- slice_sample(metadata_df, n = sample_n)
  retplot <- ggplot(metadata_df) +
    aes_string(umap1, umap2) + 
    stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE) +
    scale_fill_viridis(option = "magma") +
    coord_cartesian(expand = FALSE, xlim = c(min(metadata_df[[umap1]]), max(metadata_df[[umap1]])),
                    ylim = c(min(metadata_df[[umap2]]), max(metadata_df[[umap2]]))) 
  if (sample_n != 0) 
    retplot <- retplot + geom_point(shape = '.', col = 'white', data = data_sampled)
  return(retplot)
}

#' plot_binary_on_umap
#' 
#' Plots a binary variable on UMAP
#' 
#' @param metadata_df data frame with umap coordinates and binary variable
#' @param binary_var binary variable, must be TRUE/FALSE
#' @param grouping_strs 
#' @param umap1 string denoting UMAP_1 (optional; default: UMAP_1)
#' @param umap2 string denoting UMAP_2 (optional; default: UMAP_2)
#' @param color_str string denoting the color of detected cells (default: red)
#' @param pt_size size of points in geom_point (default: 0.5)
#' @param pt_stroke stroke of points in geom_point (default: 0.1)
#' @param pt_alpha transparency of points in geom_point (default: 0.5)
#' @param x_pos x position of count text (default: NA, sets to xmax)
#' @param y_pos y position of count text (default: NA, sets to ymax)
#' @return 
#' @export

plot_binary_on_umap <- function(metadata_df, binary_var, grouping_strs = NA, 
                                umap1 = "UMAP_1", umap2 = "UMAP_2", 
                                color_str = "red", 
                                pt_size = 0.5, pt_stroke = 0.1, pt_alpha = 0.5, 
                                x_pos = NA, y_pos = NA) {
  
  binary_var <- enquo(binary_var)
  
  # group by appropriate variables
  if (is.na(grouping_strs)) {
    grouped_df <- metadata_df %>% 
      group_by(!!binary_var) 
  
  } else if (!is.character(grouping_strs)) {
  
    stop("grouping_strs must be a character vector")
  
  } else if (length(grouping_strs == 1)) {
  
    grouped_df <- metadata_df %>% 
      group_by(!!binary_var, !!sym(grouping_strs)) 
  
  } else {
  
    grouped_df <- metadata_df %>% 
      group_by(!!binary_var, !!!syms(grouping_strs))
  
  }
  
  # calculate the number of cells
  count_df <- grouped_df %>% 
    summarise(ncells = n(), .groups = "drop")
  
  # automatically detect x_pos and y_pos
  if (is.na(x_pos)) {
    x_pos <- max(metadata_df[[umap1]])
  }
  if (is.na(y_pos)) {
    y_pos <- max(metadata_df[[umap2]])
  }
  
  p <- metadata_df %>% 
    ggplot() + 
    aes_string(umap1, umap2) + 
    aes(color = !!binary_var) + 
    geom_point(size = pt_size, stroke = pt_stroke, alpha = pt_alpha) + 
    geom_text(aes(x = x_pos, y = y_pos + 2*as.integer(!!binary_var), label = ncells), 
              data = count_df) + 
    scale_color_manual(values = c(`TRUE` = color_str, `FALSE` = "grey"))
  
  return(p)
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
