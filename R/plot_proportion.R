#' calculate_proportion
#'  
#' @param .tbl input data table
#' @param var1 first grouping variable, typically `sample_ID`
#' @param var2 second grouping variable
#' 
#' @return tibble with proportions
#' @export

calculate_proportion <- function(.tbl, var1, var2) {
  var1 <- enquo(var1)
  var2 <- enquo(var2)
  
  propr_df <- .tbl %>% 
    group_by(!!var1, !!var2) %>% 
    count() %>% 
    group_by(!!var1) %>% 
    mutate(total = sum(n), prop = n / total) %>% 
    ungroup() 
  propr_df %>% 
    pivot_wider(id_cols = !!var1, 
                names_from = !!var2, 
                values_from = prop, 
                values_fill = 0) %>% 
    pivot_longer(cols = 2:ncol(.), 
                 names_to = rlang::as_name(var2), 
                 values_to = "prop") %>% 
    left_join(unique(select(propr_df, !!var1, total)), 
              by = rlang::as_name(var1)) %>% 
    left_join(select(propr_df, !!var1, !!var2, n),
              by = c(rlang::as_name(var1),
                     rlang::as_name(var2))) %>%
    mutate(n = ifelse(is.na(n), 0, n))
}


#' plot_prop
#' 
#' Calculate the proportion and normalized proportion of entries
#' represented within a given group
#' 
#' Example: metadata %>% 
#' calculate_proportion(clustering, orig.ident) %>% 
#' plot_prop(clustering, orig.ident)
#' 
#' @param prop_df proportion data frame, usually the output of `calculate_proportion()`
#' @param var grouping var (for Seurat metadata, usually `clustering`)
#' @param norm_var var to normalize against (for Seurat metadata, usually `orig.ident`)
#' @param col_position position variable in geom_col (default: "stack"). Use ?geom_col to see other options
#' @return tibble with proportion and normalized proportion
#' @export

plot_prop <- function(prop_df, var, norm_var, col_position = "stack") {
  var <- enquo(var)
  norm_var <- enquo(norm_var)
  prop_df %>% 
    ggplot() + 
    aes(!!norm_var, prop, fill = !!var) + 
    geom_col(position = col_position) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
}

#' plot_norm_prop
#' 
#' Calculate the proportion and normalized proportion of entries
#' represented within a given group
#' 
#' Example: metadata %>% 
#' calculate_proportion(clustering, orig.ident) %>% 
#' plot_prop(clustering, orig.ident)
#' 
#' @param prop_df proportion data frame, usually the output of `calculate_proportion()`
#' @param var grouping var (for Seurat metadata, usually `clustering`)
#' @param norm_var var to normalize against (for Seurat metadata, usually `orig.ident`)
#' @param col_position position variable in geom_col (default: "stack"). Use ?geom_col to see other options
#' @return tibble with proportion and normalized proportion
#' @export
#' 
plot_norm_prop <- function(prop_df, var, norm_var, col_position = "stack") {
  var <- enquo(var)
  norm_var <- enquo(norm_var)
  prop_df %>% 
    ggplot() + 
    aes(!!var, norm_prop, fill = !!norm_var) + 
    geom_col(position = col_position) 
}
