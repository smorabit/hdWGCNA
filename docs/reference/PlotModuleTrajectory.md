# PlotModuleTrajectory

Plots module eigengene dynamics over a pseudotime trajectory.

## Usage

``` r
PlotModuleTrajectory(
  seurat_obj,
  pseudotime_col = "pseudotime",
  n_bins = 20,
  harmonized = TRUE,
  ncol = 4,
  point_size = 1,
  line_size = 1,
  se = TRUE,
  group_colors = NULL,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- pseudotime_col:

  name of the column in seurat_obj@meta.data containing pseudotime
  information. Multiple names can be passed to plot multiple
  trajectories simultaneously.

- n_bins:

  number of pseudotime bins/windows used to smooth the MEs

- harmonized:

  logical determining whether or not to use the harmonized MEs

- point_size:

  size of the points in the plot

- line_size:

  width of the lines in the plot

- se:

  logical determining whether or not to show the standard error on the
  plot

- group_colors:

  optional list of colors to differentiate multiple pseudotime
  trajectories. Must be the same length as pseudotime_col

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- n_col:

  number of columns for different modules in the plot
