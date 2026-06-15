# PlotKMEs

Plotting function to show genes by kME value in each module

## Usage

``` r
PlotKMEs(
  seurat_obj,
  n_hubs = 10,
  text_size = 2,
  ncol = 5,
  plot_widths = c(3, 2),
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- n_hubs:

  number of hub genes to display

- text_size:

  controls the size of the hub gene text

- ncol:

  number of columns to display individual plots

- plot_widths:

  the relative width between the kME rank plot and the hub gene text

- wgcna_name:

  the name of the WGCNA experiment in the seurat object
