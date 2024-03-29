% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{ModuleTFNetwork}
\alias{ModuleTFNetwork}
\title{ModuleTFNetwork}
\usage{
ModuleTFNetwork(
  seurat_obj,
  tf_name,
  tf_gene_name,
  edge.alpha = 0.75,
  cor_thresh = 0.25,
  high_color = "red",
  mid_color = "grey",
  low_color = "blue",
  slot = "data",
  size.scale = 30,
  tf_x = 0,
  tf_y = 0,
  wgcna_name = NULL
)
}
\arguments{
\item{seurat_obj}{A Seurat object}

\item{tf_name}{the Motif name for the tf}

\item{tf_gene_name}{the gene associated with this tf in the rownames(seurat_obj)}

\item{edge.alpha}{scaling factor for edge opacity in the network}

\item{cor_thresh}{threshold to plot correlation edges between modules}

\item{high_color}{color for positive correlation}

\item{mid_color}{color for zero correlation}

\item{low_color}{color for negative correlation}

\item{slot}{the slot in the seurat object to extract expression data for the tf_gene_name}

\item{size.scale}{scaling factor for the size of each node}

\item{tf_x}{x coordinate for the TF if the TF is not found in the UMAP}

\item{tf_y}{y coordinate for the TF if the TF is not foudn in the UMAP}

\item{wgcna_name}{the name of the WGCNA experiment in the seurat object}
}
\description{
Plotting the relationships between a TF and the co-expression modules
}
\examples{
ModuleTFNetwork
}
\keyword{scRNA-seq}
