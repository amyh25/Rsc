#' pb_to_deseq
#' 
#' pb to deseq
#' 
#' @param pb mtx
#' @param sample_var as string, name of variable identifying samples
#' @param vars list of strings, name of variables of relevant conditions
#' @return deseq object
#' @export

pb_to_deseq <- function(pb, sample_var, vars, 
                        formula_str = paste0("~", paste0(vars, collapse = "+"))) {
  counts_mtx <- as.matrix(pb)
  counts_metadata <- tibble(id = rownames(pb)) %>% 
    separate(id, into = c(sample_var, vars), sep = "_")
  de <- DESeqDataSetFromMatrix(
    t(counts_mtx), 
    colData = counts_metadata, 
    design = as.formula(formula_str)
  )
  de <- DESeq(de)
}

