# ModuleRegulatoryNetworkPlot

This function visualizes the regulatory network between gene modules
from a Seurat object. The plot displays transcription factor (TF)
regulatory interactions between modules, with edges representing
positive or negative regulatory links.

## Usage

``` r
ModuleRegulatoryNetworkPlot(
  seurat_obj,
  feature = "delta",
  TFs_only = TRUE,
  layout = "umap",
  umap_background = FALSE,
  max_val = 1,
  cutoff = 0,
  focus_source = NULL,
  focus_target = NULL,
  loops = TRUE,
  label_modules = TRUE,
  high_color = "orange2",
  mid_color = "white",
  low_color = "dodgerblue",
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data and WGCNA results.

- feature:

  A character string specifying the type of regulatory score to plot.
  Options are 'positive' (positive regulatory score), 'negative'
  (negative regulatory score), or 'delta' (difference between positive
  and negative scores). Default is 'delta'.

- TFs_only:

  Logical; if TRUE (default), only transcription factor (TF) genes are
  included in the plot. If FALSE, all genes are considered in the
  regulatory network.

- layout:

  A character string specifying the layout for the plot. Default is
  'umap', which arranges the modules according to UMAP coordinates.
  Other layout options from `ggraph` can also be used.

- umap_background:

  Logical; if TRUE, the UMAP of module eigengenes is plotted in the
  background (if `layout = 'umap'`). Default is FALSE.

- max_val:

  Numeric; sets the maximum absolute value for the regulatory score. Any
  values exceeding this threshold are capped. Default is 1.

- cutoff:

  Numeric; edges with absolute regulatory scores below this value are
  excluded from the plot. Default is 0.

- focus_source:

  Character vector; optionally restricts the plot to only show
  regulatory links from specified source modules. Default is NULL (all
  modules included).

- focus_target:

  Character vector; optionally restricts the plot to only show
  regulatory links targeting specified modules. Default is NULL (all
  modules included).

- loops:

  Logical; if TRUE (default), loops (self-regulatory connections) are
  shown in the plot.

- label_modules:

  Logical; if TRUE (default), module names are displayed as labels in
  the plot.

- high_color:

  Character string; the color representing high regulatory scores in the
  color gradient. Default is 'orange2'.

- mid_color:

  Character string; the color representing intermediate regulatory
  scores in the color gradient. Default is 'white'.

- low_color:

  Character string; the color representing low regulatory scores in the
  color gradient. Default is 'dodgerblue'.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- ...:

  Additional arguments passed to the
  [`ggraph::create_layout`](https://ggraph.data-imaginist.com/reference/ggraph.html)
  function.

## Value

A ggplot2 object visualizing the regulatory network, which can be
further customized or displayed using
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Details

The function visualizes the regulatory network between gene modules by
plotting transcription factor regulatory interactions. Positive and
negative regulatory scores are displayed as edges connecting nodes
(modules), with node sizes proportional to the number of genes in each
module. Users can customize the network layout (e.g., UMAP or other
graph layouts), filter edges based on score thresholds, and highlight
specific modules by using the focus parameters. Edge colors represent
the strength and direction (positive or negative) of regulatory
interactions, and loops can optionally be plotted to indicate
self-regulation.
