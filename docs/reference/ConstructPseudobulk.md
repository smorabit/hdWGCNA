# ConstructPseudobulk

Constructs a "pseudobulk" gene expression matrix summarizing the
expression levels of each gene across a grouping variable (cell types
for example) in each biological replicate.

## Usage

``` r
ConstructPseudobulk(
  seurat_obj,
  group.by,
  replicate_col,
  label_col = NULL,
  assay = "RNA",
  slot = "counts",
  layer = "counts",
  min_reps = 20,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- group.by:

  column in seurat_obj@meta.data containing grouping info, ie clusters
  or celltypes

- replicate_col:

  column in seurat_obj@meta.data denoting each replicate / sample

- label_col:

  column in seurat_obj@meta.data denoting an additional label of
  interest, for example disease status or biological sex. This is not a
  required argument and is typically only used for consensus WGCNA

- assay:

  Assay in seurat_obj containing isoform expression information.

- slot:

  Slot to extract data for aggregation. Default = 'counts'

- layer:

  Layer to extract data for aggregation. Default = 'counts'. Layer is
  used with Seurat v5 instead of slot.

- min_reps:

  The minimum number of different biological replicates allowed. Error
  will be thrown if the number of reps is too low.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

a matrix containing pseudobulk expression profiles

## Details

This function constructs pseudobulk gene expression profiles across the
provided cell groups and the provided biological replicates. We note
that low numbers of replicates are typical in single-cell and spatial
transcriptomics due to the large monetary cost of running these
experiments, and pseudobulk-ing your data for hdWGCNA is only
recommended in the case where you have a sufficient number of
replicates. Here we have set the minimum recommended number to 20. Using
fewer than 20 replicates risks the results not being reproducible or
robust, and therefore are not biologically meaningful due to spurious
noisy correlations.
