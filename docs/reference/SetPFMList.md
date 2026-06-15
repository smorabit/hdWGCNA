# SetPFMList

SetPFMList

## Usage

``` r
SetPFMList(seurat_obj, pfm_list)
```

## Arguments

- seurat_obj:

  A Seurat object

- pfm_list:

  list of pfm objects

## Examples

``` r
SetPFMList
#> function (seurat_obj, pfm_list) 
#> {
#>     if (is.null(seurat_obj@misc$motifs)) {
#>         seurat_obj@misc$motifs <- list()
#>     }
#>     seurat_obj@misc$motifs$pfm_list <- pfm_list
#>     seurat_obj
#> }
#> <bytecode: 0x55d8a3710f98>
#> <environment: namespace:hdWGCNA>
```
