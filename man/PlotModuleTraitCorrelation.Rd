% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{PlotModuleTraitCorrelation}
\alias{PlotModuleTraitCorrelation}
\title{PlotModuleTraitCorrelation}
\usage{
PlotModuleTraitCorrelation(
  seurat_obj,
  high_color = "red",
  mid_color = "grey90",
  low_color = "blue",
  label = NULL,
  label_symbol = "stars",
  plot_max = NULL,
  text_size = 2,
  text_color = "black",
  text_digits = 3,
  combine = TRUE,
  wgcna_name = NULL
)
}
\arguments{
\item{seurat_obj}{A Seurat object}

\item{plot_labels}{logical determining whether to plot the module labels#' @param wgcna_name The name of the scWGCNA experiment in the seurat_obj@misc slot}
}
\description{
Plotting function for Module Preservation statistics
}
\examples{
PlotModulePreservation
}
\keyword{scRNA-seq}