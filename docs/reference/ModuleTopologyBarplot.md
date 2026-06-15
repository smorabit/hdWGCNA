# ModuleTopologyBarplot

Plots a ranked barplot of genes in a co-expression module by
intramodular connectivity

## Usage

``` r
ModuleTopologyBarplot(
  seurat_obj,
  mod,
  features = "kME",
  plot_color = NULL,
  alpha = TRUE,
  genes_order = NULL,
  return_genes = FALSE,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- mod:

  the name of the co-expression module to plot

- features:

  specify the features to use in the barplot, 'kME' or 'degree' or
  'weighted_degree' (degree scaled to 0 or 1)

- plot_color:

  color used for the bar plot, default is the module's unique color

- alpha:

  logical indicating whether or not to add opacity to the barplot based
  on the strength (kME or degree)

- genes_order:

  a character vector of genes to plot in this specific order, this
  option will override the order_by parameter

- return_genes:

  logical indicating whether or not to return

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

ggplot object containing the ModuleTopologyBarplot

## Details

ModuleTopologyBarplot generates a barplot showing the intramodular
connectivity of each gene in a specific co-expression module. Each bar
in this plot represents a single gene, and they are ranked based on the
strength of their connections within that particular module. A custom
gene ordering can be supplied, which is helpful when comparing the
module topologies side by side with more than one dataset.
