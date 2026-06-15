# ConstructMetacells

This function takes a Seurat object and constructs averaged 'metacells'
based on neighboring cells.

## Usage

``` r
ConstructMetacells(
  seurat_obj,
  name = "agg",
  ident.group = "seurat_clusters",
  k = 25,
  reduction = "pca",
  dims = NULL,
  assay = "RNA",
  cells.use = NULL,
  slot = "counts",
  layer = "counts",
  meta = NULL,
  return_metacell = FALSE,
  mode = "average",
  max_shared = 15,
  target_metacells = 1000,
  max_iter = 5000,
  verbose = FALSE,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- name:

  A string appended to resulting metalcells. Default = 'agg'

- k:

  Number of nearest neighbors to aggregate. Default = 50

- reduction:

  A dimensionality reduction stored in the Seurat object. Default =
  'umap'

- dims:

  A vector represnting the dimensions of the reduction to use. Either
  specify the names of the dimensions or the indices. Default = NULL to
  include all dims.

- assay:

  Assay to extract data for aggregation. Default = 'RNA'

- slot:

  Slot to extract data for aggregation. Default = 'counts'. Slot is used
  with Seurat v4 instead of layer.

- layer:

  Layer to extract data for aggregation. Default = 'counts'. Layer is
  used with Seurat v5 instead of slot.

- return_metacell:

  Logical to determine if we return the metacell seurat object (TRUE),
  or add it to the misc in the original Seurat object (FALSE). Default
  to FALSE.

- mode:

  determines how to make gene expression profiles for metacells from
  their constituent single cells. Options are "average" or "sum".

- max_shared:

  the maximum number of cells to be shared across two metacells

- target_metacells:

  the maximum target number of metacells to construct

- max_iter:

  the maximum number of iterations in the metacells bootstrapping loop

- verbose:

  logical indicating whether to print additional information

- wgcna_name:

  name of the WGCNA experiment
