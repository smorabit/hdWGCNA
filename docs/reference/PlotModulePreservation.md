# PlotModulePreservation

Plotting function for Module Preservation statistics

## Usage

``` r
PlotModulePreservation(
  seurat_obj,
  name,
  statistics = "summary",
  plot_labels = TRUE,
  label_size = 4,
  mod_point_size = 4,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- name:

  The name of the module preservation analysis to plot given in
  ModulePreservation

- statistics:

  Which module preservation statistics to plot? Choices are summary,
  all, or a custom list

- plot_labels:

  logical determining whether to plot the module labels

- label_size:

  the size of the module labels

- mod_point_size:

  the size of the points in each plot

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Details

This function creates a scatter plot showing the module preservation
statistics for each module compared to the size of the module (number of
genes).
