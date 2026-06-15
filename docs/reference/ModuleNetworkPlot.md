# ModuleNetworkPlot

Visualizes the top hub genes for selected modules as a circular network
plot with one inner circle and one outer circle.

## Usage

``` r
ModuleNetworkPlot(
  seurat_obj,
  n_inner = 10,
  n_outer = 15,
  n_conns = 500,
  mods = "all",
  outdir = "ModuleNetworks",
  wgcna_name = NULL,
  plot_size = c(6, 6),
  edge.alpha = 0.25,
  edge.width = 1,
  vertex.label.cex = 1,
  vertex.size = 6,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- n_inner:

  Number of genes to plot on the inner circle

- n_outer:

  Number of genes to plot on the outer circle.

- n_conns:

  Number of gene-gene co-expression connections to plot, sorted by the
  strongest connections.

- mods:

  Names of the modules to plot. If mods = "all", all modules are
  plotted.

- outdir:

  The directory where the plots will be stored.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- plot_size:

  A vector containing the width and height of the network plots.
  example: c(5,5)

- edge.alpha:

  value between 0 and 1 determining the alpha (transparency) scaling
  factor for the network edges

- edge.width:

  value determining the width of the network edges

- vertex.label.cex:

  vertex label font size

- vertex.size:

  vertex size
