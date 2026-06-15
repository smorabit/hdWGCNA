# MetaspotsByGroups

Computes metaspots in a given Seurat object containing spatial
transcriptomics data.

## Usage

``` r
MetaspotsByGroups(
  seurat_obj,
  group.by = c("seurat_clusters"),
  ident.group = "seurat_clusters",
  assay = "Spatial",
  slot = "counts",
  layer = "counts",
  mode = "sum",
  min_spots = 50,
  wgcna_name = NULL
)
```

## Arguments

- group.by:

  A character vector of Seurat metadata column names representing groups
  for which metacells will be computed.

- assay:

  Assay to extract data for aggregation. Default = 'Spatial'

- slot:

  Slot to extract data for aggregation. Default = 'counts'

- layer:

  Layer to extract data for aggregation. Default = 'counts'. Layer is
  used with Seurat v5 instead of slot.

- mode:

  determines how to make gene expression profiles for metacells from
  their constituent single cells. Options are "average" or "sum".

- min_spots:

  the minimum number of spots in a particular grouping to construct
  metaspots

- wgcna_name:

  name of the WGCNA experiment

- cur_seurat:

  A Seurat object
