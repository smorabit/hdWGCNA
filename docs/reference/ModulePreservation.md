# ModulePreservation

Computes module preservation statistics in a query dataset for a given
reference dataset

## Usage

``` r
ModulePreservation(
  seurat_obj,
  seurat_ref,
  name,
  n_permutations = 500,
  parallel = FALSE,
  seed = 12345,
  gene_mapping = NULL,
  genome_ref_col = NULL,
  genome_query_col = NULL,
  return_raw = FALSE,
  wgcna_name = NULL,
  wgcna_name_ref = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- seurat_ref:

  A Seurat object serving as the reference for the module preservation
  analysis

- name:

  The name to give the module preservation analysis.

- n_permutations:

  Number of permutations for the module preservation test.

- parallel:

  logical determining whether to run preservation analysis in parallel

- seed:

  random seed for the permutation analysis.

- gene_mapping:

  a dataframe containing gene name mappings between the query and the
  reference dataset. One column should have the gene name in the query
  dataset, and anotehr column should have the corresponding gene name in
  the reference dataset.

- genome_ref_col:

  the column name containing the gene names for the reference dataset

- genome_query_col:

  the column name containing the gene names for the query dataset

- return_raw:

  if TRUE, returns the module preservation statistics, else returns
  seurat_obj with the stats added to the hdWGCNA experiment.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- wgcna_name_ref:

  The name of the hdWGCNA experiment in the seurat_ref@misc slot

## Details

ModulePreservation performs a statistical test to assess the
preservation of co-expression modules identified in one dataset in an
independent dataset. This method is originally described by Langfelder
et al. in the 2011 paper "Is My Network Module Preserved and
Reproducible?". This method can be used to assess biological differences
in networks, as well as technical differences / reproducibility across
different batches.
