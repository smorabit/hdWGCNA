# GetMELoadings

Function to retrieve module eigengens from Seurat object.

## Usage

``` r
GetMELoadings(seurat_obj, harmonized = TRUE, wgcna_name = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object

- harmonized:

  logical indicating whether MEs have been harmonized

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot
