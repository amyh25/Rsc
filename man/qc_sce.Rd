% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk_create.R
\name{qc_sce}
\alias{qc_sce}
\title{qc_sce}
\usage{
qc_sce(sce_obj, sample_var, threshold = 10)
}
\arguments{
\item{sce_obj}{SingleCellExperiment object}

\item{sample_var}{as string, name of variable identifying samples}

\item{threshold}{min counts per gene threshold (default: 10)}
}
\value{
sce object
}
\description{
QC for sce
}
