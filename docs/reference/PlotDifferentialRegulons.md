# PlotDifferentialRegulons

Function to visualize differential TF regulon activity between two
groups of cells based on regulon scores. The plot shows the average log
fold-change (logFC) for positive and negative regulon scores in a
scatter plot, with the points colored by their associated module and
optionally labeled with the most differentially active regulons.

## Usage

``` r
PlotDifferentialRegulons(
  seurat_obj,
  dregs,
  n_label = 10,
  logfc_thresh = 0.1,
  lm = TRUE,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object that contains the WGCNA experiment and TF regulon
  scores.

- dregs:

  Dataframe output from the `FindDifferentialRegulons` function
  containing differential regulon results.

- n_label:

  Integer specifying the number of top up- and down-regulated regulons
  to label on the plot. If 'all', all significant regulons will be
  labeled. Default = 10.

- logfc_thresh:

  Numeric threshold for labeling and annotating regulons based on logFC
  values. Default = 0.1.

- lm:

  Logical indicating whether to plot a linear regression line for the
  logFC values of positive vs. negative regulons. Default = TRUE.

- wgcna_name:

  The name of the WGCNA experiment in the `seurat_obj@misc` slot.
  Default is the active WGCNA experiment.

## Value

A ggplot object representing the differential regulon scatter plot.

## Details

This function generates a scatter plot where each point represents a TF
regulon, with the x-axis showing the average log fold-change for
positive regulon scores and the y-axis showing the average log
fold-change for negative regulon scores. Points are colored by their
module assignment, and the size of the points reflects the module
eigengene correlation (kME) value of the TF.

Differentially expressed regulons can be highlighted, and the most
significantly up- and down-regulated regulons can be labeled. Additional
options allow adding a linear regression line and controlling label
density based on logFC thresholds.
