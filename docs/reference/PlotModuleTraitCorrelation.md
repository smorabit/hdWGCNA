# PlotModuleTraitCorrelation

Plotting function for Module Preservation statistics

## Usage

``` r
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
```

## Arguments

- seurat_obj:

  A Seurat object

- high_color:

  color for positive correlation

- mid_color:

  color for zero correlation

- low_color:

  color for negative correlation

- label:

  logical determining whether to add p-val label in each cell of the
  heatmap

- label_symbol:

  show the labels as 'stars' or as 'numeric'

- plot_max:

  maximum value of correlation to show on the colorbar

- text_size:

  size of the labels

- text_color:

  color of the text labels

- text_digits:

  how many digits to show in the text labels

- combine:

  logical determining whether to plot as one combined plot (TRUE) or to
  return individual plots as a list (FALSE)
