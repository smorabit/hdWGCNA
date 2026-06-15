# ModulePreservationNetRep

Computes module preservation statistics in a query dataset for a given
reference dataset

## Usage

``` r
ModulePreservationNetRep(
  seurat_query,
  seurat_ref,
  name,
  n_permutations = 10000,
  n_threads = 8,
  gene_mapping = NULL,
  genome_ref_col = NULL,
  genome_query_col = NULL,
  TOM_use = NULL,
  ref_power = 1,
  query_power = 1,
  wgcna_name = NULL,
  wgcna_name_ref = NULL,
  ...
)
```

## Arguments

- seurat_query:

  A Seurat object

- seurat_ref:

  A Seurat object serving as the reference for the module preservation
  analysis

- name:

  The name to give the module preservation analysis.

- n_permutations:

  Number of permutations for the module preservation test.

- n_threads:

  Number of parallel threads to use for the module preservation test

- gene_mapping:

  a dataframe containing gene name mappings between the query and the
  reference dataset. One column should have the gene name in the query
  dataset, and anotehr column should have the corresponding gene name in
  the reference dataset.

- genome_ref_col:

  the column name containing the gene names for the reference dataset

- genome_query_col:

  the column name containing the gene names for the query dataset

- TOM_use:

  The name of the hdWGCNA experiment containing the TOM that will be
  used for plotting

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- wgcna_name_ref:

  The name of the hdWGCNA experiment in the seurat_ref@misc slot

## Value

A Seurat Object containing the module preservation statistics from
NetRep.

## Details

ModulePreservationNetRep performs a statistical test to assess the
preservation of co-expression modules identified in one dataset in an
independent dataset using the NetRep implementation rather than the
WGCNA implementation. This method is originally described by Langfelder
et al. in the 2011 paper "Is My Network Module Preserved and
Reproducible?". The NetRep implementation of the module preservation
test is described by Ritchie et al. in the 2016 paper "A Scalable
Permutation Approach Reveals Replication and Preservation Patterns of
Networks Modules in Large Datasets" This method can be used to assess
biological differences in networks, as well as technical differences /
reproducibility across different batches. Note that the outputs of
ModulePreservation and ModulePreservationNetRep are not identical. This
function requires separate installation of the NetRep package, which is
not directly included as a dependency in hdWGCNA.
