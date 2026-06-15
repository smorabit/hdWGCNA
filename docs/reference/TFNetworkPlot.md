# TFNetworkPlot

Plots a network of TFs and predicted target genes.

## Usage

``` r
TFNetworkPlot(
  seurat_obj,
  selected_tfs,
  depth = 2,
  edge_weight = "Cor",
  cutoff = 0.01,
  color_cutoff = 0.75,
  target_type = "both",
  use_regulons = TRUE,
  label_genes = NULL,
  label_TFs = 1,
  no_labels = FALSE,
  TFs_only = FALSE,
  high_color = "orange2",
  mid_color = "white",
  low_color = "dodgerblue",
  node_colors = c("black", "darkorchid4", "mediumpurple2"),
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- selected_tfs:

  A list of TFs

- depth:

  Number of layers to extend the TF network from the selected_tfs. For
  example, if depth=2 (default), the target genes of selected_tfs are
  shown, and additional target genes are shown for other TFs that are
  target genes of the original selected_tfs.

- edge_weight:

  Attribute to use to color the network edges. "Cor" to show the pearson
  correlation, or "Gain" to show the importance score from the XGBoost
  model.

- cutoff:

  Cutoff for the edge weights, links below this value will not be
  included.

- color_cutoff:

  Maximum value for the colorscale for the edge weights.

- target_type:

  Type of target genes to show in the network. "Positive" shows genes
  with expression positively correlated with selected_tfs, "Negative"
  shows genes that are negatively correlated, and "Both" shows both
  (default).

- use_regulons:

  Logical indicating whether to use the TF-gene links defined in the TF
  Regulons (TRUE, default), or all putative TF-gene links (FALSE).

- label_genes:

  List of target genes to label in the network plot.

- label_TFs:

  The depth level in the network to plot the names of TFs.

- no_labels:

  Logical indicating whether or not to remove all labels.

- TFs_only:

  Logical indicating whether or not to use only TFs in the plot or to
  use TFs and other genes.

- high_color:

  Color corresponding to positive edge weights

- mid_color:

  Color corresponding to edge weights close to 0

- low_color:

  Color corresponding to negative edge weights

- node_colors:

  List of color names for the TFs and genes at each depth. For example,
  if depth=2, you should supply a list of 3 colors corresponding to the
  selected_tfs (depth 0), the primary target genes (depth 1), and the
  secondary target genes (depth 2).

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

ggplot object containing the TFNetworkPlot

## Details

TFNetworkPlot plots a network of TFs and their predicted target genes,
based on the results from ConstructTFNetwork.
