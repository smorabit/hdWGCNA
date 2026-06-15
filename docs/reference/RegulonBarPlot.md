# RegulonBarPlot

Plots the top target genes within a specific TF regulon

## Usage

``` r
RegulonBarPlot(
  seurat_obj,
  selected_tf,
  cutoff = 0.2,
  top_n = Inf,
  TFs_only = FALSE,
  high_color = "orange2",
  mid_color = "white",
  low_color = "dodgerblue",
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- cutoff:

  Cutoff for the edge weights, links below this value will not be
  included.

- top_n:

  Number of top and bottom target genes to include in the plot.

- TFs_only:

  Logical indicating whether or not to use only TFs in the plot or to
  use TFs and other genes.

- high_color:

  Color corresponding to positive edge weights

- mid_color:

  Color corresponding to edge weights close to 0

- low_color:

  Color corresponding to negative edge weights

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- selected_tfs:

  A TF whose target genes will be shown

## Value

ggplot object containing the RegulonBarPlot

## Details

RegulonBarPlot creates a bar plot showing the top target genes for a
particular TF based on the TF regulatory network analysis. Target genes
are ranked by the predicted interaction strength from XGBoost (Gain)
multiplied by the sign of the correlation between the TF and the target
gene.
