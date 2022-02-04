% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_proportion.R
\name{plot_prop}
\alias{plot_prop}
\title{plot_prop}
\usage{
plot_prop(prop_df, var, norm_var)
}
\arguments{
\item{prop_df}{proportion data frame, usually the output of \code{calculate_proportion()}}

\item{var}{grouping var (for Seurat metadata, usually \code{clustering})}

\item{norm_var}{var to normalize against (for Seurat metadata, usually \code{orig.ident})}
}
\value{
tibble with proportion and normalized proportion
}
\description{
Calculate the proportion and normalized proportion of entries
represented within a given group
}
\details{
Example: metadata \%>\%
calculate_proportion(clustering, orig.ident) \%>\%
plot_prop(clustering, orig.ident)
}