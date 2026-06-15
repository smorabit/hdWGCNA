# ResetModuleNames

Reset the uname of each hdWGCNA module

## Usage

``` r
ResetModuleNames(
  seurat_obj,
  new_name = "M",
  reset_levels = FALSE,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- new_name:

  string containing the base name to re-name the modules

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot
