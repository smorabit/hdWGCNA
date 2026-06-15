# FindMajorIsoforms

Finds the set of "major isoforms" for each gene in each cell population.
The set of major isoforms consists of the isoforms accounting for a
desired proportion of a gene's expression. Pseudobulk replicates in each
sample and cell population are used for this calculation.

## Usage

``` r
FindMajorIsoforms(
  seurat_obj,
  group.by,
  replicate_col,
  isoform_delim = "[.]",
  proportion_thresh = 0.8,
  low_thresh = 25,
  assay = "iso",
  slot = "counts",
  cluster_markers = NULL,
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

- isoform_delim:

- proportion_thresh:

  desired proportion of expression to define the set of major isoforms.
  Default = 0.8.

- low_thresh:

  lower bound of expression level for considering an isoform as part of
  the major isoform set.

- assay:

  Assay in seurat_obj containing isoform expression information.

- slot:

  Slot in seurat_obj, default to counts slot.

- cluster_markers:

  Cell population marker gene table from Seurat FindAllMarkers for the
  same cell populations specified in group.by. Optional parameter, will
  exclude isoforms that are not from marker genes.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

a list of major isoforms for each cell population

## Details

FindMajorIsoforms computes the set of major isoforms in a given Seurat
object that contains isoform-level expression information. First,
pseudobulk replicates are computed for the given cell populations and
samples present in the Seurat object. For each gene in each cell
population, we rank the gene's isoforms by expression level and take the
top expressing isoforms that make up the desired proportion of the
gene's total expression, making sure to exclude any very lowly expressed
isoforms.

Optionally, the user may supply a marker gene table for each cell
population (formatted like the output of Seurat FindAllMarkers), and
then the algorithm will only return major isoforms of the given marker
genes.
