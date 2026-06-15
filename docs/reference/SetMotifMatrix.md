# SetMotifMatrix

SetMotifMatrix

## Usage

``` r
SetMotifMatrix(seurat_obj, tf_match)
```

## Arguments

- seurat_obj:

  A Seurat object

- tf_match:

  matrix containing tf-promoter matches

## Examples

``` r
SetMotifMatrix
#> function (seurat_obj, tf_match) 
#> {
#>     if (is.null(seurat_obj@misc$motifs)) {
#>         seurat_obj@misc$motifs <- list()
#>     }
#>     seurat_obj@misc$motifs$tf_match_matrix <- tf_match
#>     seurat_obj
#> }
#> <bytecode: 0x55d8a2b247c8>
#> <environment: namespace:hdWGCNA>
```
