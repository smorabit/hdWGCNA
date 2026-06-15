# SetMotifs

SetMotifs

## Usage

``` r
SetMotifs(seurat_obj, motif_df)
```

## Arguments

- seurat_obj:

  A Seurat object

- motif_df:

  dataframe containing info about the motifs being analyzed

## Examples

``` r
SetMotifs
#> function (seurat_obj, motif_df) 
#> {
#>     if (is.null(seurat_obj@misc$motifs)) {
#>         seurat_obj@misc$motifs <- list()
#>     }
#>     seurat_obj@misc$motifs$motif_df <- motif_df
#>     seurat_obj
#> }
#> <bytecode: 0x55d89a313ab8>
#> <environment: namespace:hdWGCNA>
```
