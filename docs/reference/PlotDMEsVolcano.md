# PlotDMEsVolcano

Plotting function for the results of FindDMEs and FindAllDMEs.

## Usage

``` r
PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  plot_labels = TRUE,
  mod_point_size = 4,
  label_size = 4,
  show_cutoff = TRUE,
  wgcna_name = NULL,
  xlim_range = NULL,
  ylim_range = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object containing the WGCNA analysis in the @misc slot.

- DMEs:

  A dataframe output from FindDMEs or FindAllDMEs containing DME
  results.

- plot_labels:

  Logical, determines whether to plot the module labels on the volcano
  plot. Default is TRUE.

- mod_point_size:

  Numeric, the size of the points on the volcano plot. Default is 4.

- label_size:

  Numeric, the size of the module labels on the plot. Default is 4.

- show_cutoff:

  Logical, determines whether to plot the significance cutoff. Set this
  to FALSE if using facet_wrap. Default is TRUE.

- wgcna_name:

  Character, the name of the hdWGCNA experiment in the `seurat_obj@misc`
  slot. Default is NULL, in which case it pulls the active WGCNA
  experiment from `seurat_obj@misc$active_wgcna`.

- xlim_range:

  A numeric vector of length 2 specifying the x-axis limits for the log2
  fold change. Default is NULL, which automatically calculates limits
  based on the data.

- ylim_range:

  A numeric vector of length 2 specifying the y-axis limits for the
  -log10(p-value). Default is NULL, which automatically calculates
  limits based on the data.

## Value

A ggplot object containing the volcano plot for the DME results.

## Details

This function generates a volcano plot for differential module
expression (DME) analysis results. It can handle both two-group
comparisons (using the output of `FindDMEs`) and one-vs-all comparisons
(using the output of `FindAllDMEs`).

## Examples

``` r
# Example usage:
# Assuming `seurat_obj` is your Seurat object and `DMEs` is the output from FindDMEs
PlotDMEsVolcano(seurat_obj, DMEs, wgcna_name = "MG")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'na.omit': object 'DMEs' not found
```
