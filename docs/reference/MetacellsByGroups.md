# MetacellsByGroups

This function takes a Seurat object and constructs averaged 'metacells'
based on neighboring cells in provided groupings, such as cluster or
cell type.

## Usage

``` r
MetacellsByGroups(
  seurat_obj,
  group.by = c("seurat_clusters"),
  ident.group = "seurat_clusters",
  k = 25,
  reduction = "pca",
  dims = NULL,
  assay = NULL,
  slot = "counts",
  layer = "counts",
  mode = "average",
  cells.use = NULL,
  min_cells = 100,
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

- group.by:

  A character vector of Seurat metadata column names representing groups
  for which metacells will be computed.

- k:

  Number of nearest neighbors to aggregate. Default = 50

- reduction:

  A dimensionality reduction stored in the Seurat object. Default =
  'pca'

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

- mode:

  determines how to make gene expression profiles for metacells from
  their constituent single cells. Options are "average" or "sum".

- min_cells:

  the minimum number of cells in a particular grouping to construct
  metacells

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

- name:

  A string appended to resulting metalcells. Default = 'agg'

## Value

seurat_obj with a metacell seurat_obj stored in the specified WGCNA
experiment

## Details

MetacellsByGroups merges transcriptomically similar cells into
"metacells". Given a dimensionally-reduced representation of the input
dataset, this algorithm first uses KNN to identify similar cells. A
bootstrapped sampling procedure is then used to group together similar
cells until convergence is reached. Importantly, this procedure is done
in a context-specific manner based on the provided group.by parameters.
Typically this means that metacells will be constructed separately for
each biological replicate, cell type or cell state, disease condition,
etc. The metacell representation is considerably less sparse than the
original single-cell dataset, which is preferable for co-expression
network analysis or othter analyses that rely on correlations.
