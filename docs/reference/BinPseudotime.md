# BinPseudotime

Makes evenly-spaced bins of cells along a pseudotime trajectory.

## Usage

``` r
BinPseudotime(seurat_obj, n_bins = 50, pseudotime_col = "pseudotime")
```

## Arguments

- seurat_obj:

  A Seurat object

- n_bins:

  number of hub genes to use in the UMAP computation

- pseudotime_col:

  name of the column in seurat_obj@meta.data containing pseudotime
  information
