# FindDifferentialRegulons

Function to compare TF regulon scores between two sets of cell barcodes.

## Usage

``` r
FindDifferentialRegulons(
  seurat_obj,
  barcodes1,
  barcodes2,
  assay = "RNA",
  slot = "data",
  layer = "data",
  wgcna_name = NULL,
  test.use = "wilcox",
  only.pos = FALSE,
  logfc.threshold = 0,
  min.pct = 0,
  verbose = FALSE,
  pseudocount.use = 0,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- barcodes1:

  character vector containing cell barcodes for the first group to test.
  Positive fold-change means up-regulated in this group.

- barcodes2:

  character vector containing cell barcodes for the second group to
  test. Negative fold-change means up-regulated in this group.

- assay:

  Assay to extract data for aggregation. Default = 'RNA'

- slot:

  Slot to extract data for aggregation. Default = 'counts'. Slot is used
  with Seurat v4 instead of layer.

- layer:

  Layer to extract data for aggregation. Default = 'counts'. Layer is
  used with Seurat v5 instead of slot.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- ...:

  Additional parameters for the Seurat FindMarkers function

## Value

A dataframe contaning differential regulon results

## Details

FindDifferentialRegulons compares two groups based on their TF regulon
scores. Three comparisons are made for different sets of features:
positive regulon scores, negative regulon scores, and gene expression.
The same settings for the test will be used for all three tests.
