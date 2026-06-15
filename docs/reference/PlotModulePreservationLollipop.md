# PlotModulePreservationLollipop

Plots a ranked lollipop plot of co-expression modules based on the
results of module preservation analysis.

## Usage

``` r
PlotModulePreservationLollipop(
  seurat_obj,
  name,
  features = NULL,
  fdr = TRUE,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- name:

  The name to give the module preservation analysis.

- features:

  The name of the module preservation features to plot.

- fdr:

  logical indicating whether or not to plot FDR-corrected p-values

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

ggplot object containing the PlotModulePreservationLollipop

## Details

PlotModulePreservationLollipop generates a lollipop plot showing module
preservation results. If the module preservation test was performed
using the WGCNA method, the statistic that will be shown is the
Z-summary preservation statistic. If the analysis was performed using
NetRep, then the statistic that will be shown is the FDR corrected
averaged p-values from the module preservation permutation test.
