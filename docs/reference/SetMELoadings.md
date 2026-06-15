# SetMELoadings

SetMELoadings

## Usage

``` r
SetMELoadings(seurat_obj, loadings, harmonized = TRUE, wgcna_name = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object

- loadings:

  named numeric vector with eigengene loadings

- harmonized:

  logical indicating whether MEs have been harmonized

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot
