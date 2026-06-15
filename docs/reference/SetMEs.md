# SetMEs

SetMEs

## Usage

``` r
SetMEs(seurat_obj, MEs, harmonized = TRUE, wgcna_name = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object

- MEs:

  dataframe or matrix containing module eigengenes

- harmonized:

  logical indicating whether MEs have been harmonized

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot
