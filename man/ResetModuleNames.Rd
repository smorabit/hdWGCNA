% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters_and_setters.R
\name{ResetModuleNames}
\alias{ResetModuleNames}
\title{ResetModuleNames}
\usage{
ResetModuleNames(
  seurat_obj,
  new_name = "M",
  reset_levels = FALSE,
  wgcna_name = NULL
)
}
\arguments{
\item{seurat_obj}{A Seurat object}

\item{new_name}{string containing the base name to re-name the modules}

\item{wgcna_name}{The name of the hdWGCNA experiment in the seurat_obj@misc slot}
}
\description{
Reset the uname of each hdWGCNA module
}
\keyword{scRNA-seq}
