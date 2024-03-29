% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_volcano_skeleton}
\alias{plot_volcano_skeleton}
\title{plot_volcano_skeleton}
\usage{
plot_volcano_skeleton(
  .data,
  fc,
  pval,
  remove_inf_fc = TRUE,
  pt.size = 0.5,
  fc_thresh = 1
)
}
\arguments{
\item{.data}{data frame with p values and fold change}

\item{fc}{name of fold change variable}

\item{pval}{name of p-val variable}

\item{remove_inf_fc}{(optional; default: TRUE) remove zero and infinite fold change values}

\item{pt.size}{(optional; default: 0.5) size of points}

\item{fc_thresh}{(optional; default: 1) fold change threshold to draw lines at}
}
\value{
ggplot of volcano plot
}
\description{
Plots skeleton of volcano plot
}
