% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_gene_str}
\alias{plot_gene_str}
\title{plot_gene_str}
\usage{
plot_gene_str(
  gene_str,
  genes_df,
  umap_coords,
  join_col_str = "cell",
  UMAP_1 = "UMAP_1",
  UMAP_2 = "UMAP_2",
  stroke_size = 0.1,
  pt_size = 0.5
)
}
\arguments{
\item{gene_str}{gene name, as string}

\item{genes_df}{data frame with at least two columns, one with gene values and one with obs IDs to join to umap_coords}

\item{umap_coords}{data frame with at least three columns, two with UMAP coords and one with obs IDs to join to genes_df}

\item{UMAP_1}{string denoting UMAP_1 (optional; default: UMAP_1)}

\item{UMAP_2}{string denoting UMAP_2 (optional; default: UMAP_2)}

\item{stroke_size}{stroke size (optional; default: 0.1)}

\item{pt_size}{point size (optional; default: 0.5)}
}
\value{
a ggplot object
}
\description{
Plots genes, with points ordered by expression of gene value
}
