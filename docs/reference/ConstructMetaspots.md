# ConstructMetaspots

Computes metaspots in a given Seurat object containing spatial
transcriptomics data. This function is called by MetaspotsByGroups and
should NOT be run directly!

## Usage

``` r
ConstructMetaspots(
  cur_seurat,
  mode = "sum",
  assay = "Spatial",
  slot = "counts",
  layer = "counts"
)
```

## Arguments

- cur_seurat:

  A Seurat object

- mode:

  "sum" or "average"

- assay:

  Assay to extract data for aggregation. Default = 'Spatial'

- slot:

  Slot to extract data for aggregation. Default = 'counts'

- layer:

  Layer to extract data for aggregation. Default = 'counts'. Layer is
  used with Seurat v5 instead of slot.
