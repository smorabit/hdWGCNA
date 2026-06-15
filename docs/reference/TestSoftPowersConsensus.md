# TestSoftPowersConsensus

Compute the scale-free topology model fit for different soft power
thresholds separately for each input dataset

## Usage

``` r
TestSoftPowersConsensus(
  seurat_obj,
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  use_metacells = TRUE,
  networkType = "signed",
  corFnc = "bicor",
  setDatExpr = FALSE,
  group.by = NULL,
  group_name = NULL,
  multi.group.by = NULL,
  multi_groups = NULL,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- powers:

  numeric vector specifying soft powers to test

- use_metacells:

  logical flag for whether to use the metacell expression matrix

- networkType:

  The type of network to use for network analysis. Options are "signed"
  (default), "unsigned", or "signed hybrid". This should be consistent
  with the network chosen for ConstructNetwork

- corFnc:

  Correlation function for the gene-gene correlation adjacency matrix.

- setDatExpr:

  logical flag indicating whether to run setDatExpr.

- group.by:

  A string containing the name of a column in the Seurat object with
  cell groups (clusters, cell types, etc). If NULL (default), hdWGCNA
  uses the Seurat Idents as the group.

- group_name:

  A string containing a group present in the provided group.by column or
  in the Seurat Idents. A character vector can be provided to select
  multiple groups at a time.

- multi.group.by:

  A string containing the name of a column in the Seurat object with
  groups for consensus WGCNA (dataset, sample, condition, etc)

- multi_groups:

  A character vecrtor containing the names of groups to select

- ...:

  additional parameters passed to SetDatExpr

## Examples

``` r
# TestSoftPowers(pbmc)
```
