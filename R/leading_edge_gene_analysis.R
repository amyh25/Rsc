#' get_leading_edge_df
#' 
#' Creates a "leading edge df" for input into other leading edge
#' analysis functions
#' 
#' @param rank_vec named vector of values
#' @param gs_vec vector of genes in gene set
#' @return leading edge df, columns: name, value, in_gs, rank, running_es, direction, is_leading_edge
#' @export

get_leading_edge_df <- function(rank_vec, gs_vec) {
  
  tibble(name = names(rank_vec), value = rank_vec) %>% 
    arrange(desc(value)) %>% 
    mutate(in_gs = name %in% gs_vec) %>% 
    mutate(rank = row_number()) %>% 
    mutate(score = ifelse(in_gs, 1 / sum(in_gs), -1 / (max(rank) - sum(in_gs)))) %>% 
    mutate(running_es = cumsum(score)) %>% 
    mutate(max_es = max(running_es), min_es = min(running_es), 
           max_abs_es = max(abs(running_es))) %>% 
    mutate(direction = ifelse(max_abs_es == max_es, "up", "dn")) %>% 
    mutate(inflection_rank = case_when(
      direction == "up" & running_es == max_es ~ rank, 
      direction == "dn" & running_es == min_es ~ rank, 
      TRUE ~ NA_integer_
    )) %>% 
    tidyr::fill(inflection_rank, .direction = "downup") %>%
    mutate(is_leading_edge = case_when(
      direction == "up" ~ rank <= inflection_rank, 
      direction == "dn" ~ rank >= inflection_rank
    )) %>% 
    select(name, value, in_gs, rank, running_es, direction, is_leading_edge)
}

#' plot_gsea_curve
#' 
#' Given a "leading edge df" and title, 
#' plots GSEA curve with power = 0
#' 
#' @param leading_edge_df dataframe output by get_leading_edge_df
#' @param title string; title of gene set
#' @param strwidth (optional; default 30) width of title strwrap
#' @param save (optional) TRUE to save the plot to the output directory
#' @param output_dir (optional) output directory if saving
#' @return ggplot
#' @export

plot_gsea_curve <- function(leading_edge_df, title, 
                            strwidth = 30, save = FALSE, output_dir = "") {
  
  pretty_title <- str_replace_all(title, "_", " ") %>% 
    str_wrap(width = strwidth)
  
  n_in_gs <- nrow(leading_edge_df %>% filter(in_gs))
  n_leg <- nrow(leading_edge_df %>% filter(in_gs) %>% filter(is_leading_edge))
  
  p <- leading_edge_df %>% 
    filter(in_gs) %>% 
    ggplot() + 
    aes(rank, running_es, color = is_leading_edge) + 
    geom_hline(yintercept = 0) + 
    geom_point(show.legend = FALSE) + 
    scale_color_manual(values = boolean_palette) +
    labs(title = pretty_title, 
         subtitle = paste0(n_in_gs, " genes in gene set; ", 
                           n_leg, " leading edge genes") %>% 
           str_wrap(width = strwidth), 
         y = "Running Enrichment Score", 
         x = "Gene Rank") + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  if (save) {
    ggsave(file.path(output_dir, paste0("gsea.curve_", title, ".pdf")), 
           height = 5, width = 5*1.618, 
           plot = p)
  }
  
  return(p)
}

#' plot_top_leg
#' 
#' Given a "leading edge df" and title, 
#' plots a bar plot with the top leading edge genes
#' 
#' @param leading_edge_df dataframe output by get_leading_edge_df
#' @param title string; title of gene set
#' @param topn (optional; default 20) the number of top leading edge genes to plot
#' @param save (optional) TRUE to save the plot to the output directory
#' @param output_dir (optional) output directory if saving
#' @return ggplot
#' @export

plot_top_leg <- function(leading_edge_gene_df, title, 
                         topn = 20, save = FALSE, output_dir = "") {
  
  if (unique(leading_edge_gene_df$direction) == "up") {
    plot_data <- leading_edge_gene_df %>% 
      filter(in_gs & is_leading_edge) %>% 
      slice_min(rank, n = topn) %>% 
      mutate(name = fct_reorder(name, -value)) 
  } else if (unique(leading_edge_gene_df$direction) == "dn") {
    plot_data <- leading_edge_gene_df %>% 
      filter(in_gs & is_leading_edge) %>% 
      slice_max(rank, n = topn) %>% 
      mutate(name = fct_reorder(name, -value)) 
  } else {
    stop("Leading edge gene data frame missing direction")
  }
  
  p <- plot_data %>% 
    ggplot() + 
    aes(name, value, fill = is_leading_edge) + 
    geom_col(show.legend = FALSE) + 
    coord_flip() + 
    scale_fill_manual(values = boolean_palette) + 
    ylab("-log10(p-val)") + 
    xlab("gene") + 
    ggtitle(paste0("Top ", topn, " Leading Edge Genes"))
  
  if (save) {
    ggsave(file.path(output_dir, paste0("gsea.top.leg_", title, ".pdf")), 
           height = 5, width = 5, 
           plot = p)
  }
  
  return(p)
  
}

#' plot_gsea_leg_plots
#' 
#' Given a "leading edge df" and title, 
#' plots GSEA curve with power = 0 and
#' plots a bar plot with the top leading edge genes
#' using `plot_grid`
#' 
#' @param leading_edge_df dataframe output by get_leading_edge_df
#' @param title string; title of gene set
#' @param topn (optional; default 20) the number of top leading edge genes to plot
#' @param strwidth (optional; default 30) width of title strwrap
#' @return grob of ggplots
#' @export

plot_gsea_leg_plots <- function(leading_edge_gene_df, title = "", 
                                strwidth = 30, topn = 20) {
  plot_grid(
    plot_gsea_curve(leading_edge_gene_df, strwidth = strwidth, title = title), 
    plot_top_leg(leading_edge_gene_df, topn = topn, title = title)
  )
}