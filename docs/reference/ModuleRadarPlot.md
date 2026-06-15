# ModuleRadarPlot

Plots the expression level (module eigengene) of each co-expression
module for different groups as a radar plot.

## Usage

``` r
ModuleRadarPlot(
  seurat_obj,
  group.by = NULL,
  barcodes = NULL,
  combine = TRUE,
  ncol = 4,
  features = "hMEs",
  wgcna_name = NULL,
  fill = TRUE,
  draw.points = FALSE,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- group.by:

  the column name of the selected comparison in the DMEs dataframe

- barcodes:

  A list of barcodes from colnames(seurat_obj) which will be used to
  subset the data before plotting.

- combine:

  logical indicating whether or not to combine plots using patchwork

- ncol:

  The number of columns for the combined plot if patchwork is being used

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- ...:

  additional parameters for ggradar

## Value

ggplot object containing the ModuleRadarPlot

## Details

ModuleRadarPlot visualizes the expression level (module eigengene, ME)
of each co-expression module for different groups in a radial coordinate
system. The ME for each module is averaged for each group, scaled, and
then is plotted radially. The resulting plots help us to interpret which
cell groups (clusters, cell types, etc) are expressing each module.
