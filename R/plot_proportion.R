#' calculate_proportion
#' 
#' Calculate the proportion and normalized proportion of entries
#' represented within a given group
#' 
#' Example: metadata %>% 
#' calculate_proportion(clustering, orig.ident) %>% 
#' plot_prop(clustering, orig.ident)
#' 
#' @param .data input data frame, where each row is an observation
#' @param var grouping var (for Seurat metadata, usually `clustering`)
#' @param norm_var var to normalize against (for Seurat metadata, usually `orig.ident`)
#' @return tibble with proportion and normalized proportion
#' @export

calculate_proportion <- function(.data, var, norm_var) {
  var <- enquo(var)
  norm_var <- enquo(norm_var)
  .data %>% 
    group_by(!!var, !!norm_var) %>% 
    summarise(count = n()) %>% 
    group_by(!!norm_var) %>% 
    mutate(total = sum(count), prop = count / total) %>% 
    group_by(!!var) %>% 
    mutate(total_prop = sum(prop), norm_prop = prop / total_prop) 
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
#' @return tibble with proportion and normalized proportion
#' @export

plot_prop <- function(prop_df, var, norm_var) {
  var <- enquo(var)
  norm_var <- enquo(norm_var)
  prop_df %>% 
    ggplot() + 
    aes(!!norm_var, prop, fill = !!var) + 
    geom_col() + 
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
#' @return tibble with proportion and normalized proportion
#' @export
#' 
plot_norm_prop <- function(prop_df, var, norm_var) {
  var <- enquo(var)
  norm_var <- enquo(norm_var)
  prop_df %>% 
    ggplot() + 
    aes(!!var, norm_prop, fill = !!norm_var) + 
    geom_col() 
}
