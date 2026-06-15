# SetMotifTargets

SetMotifTargets

## Usage

``` r
SetMotifTargets(seurat_obj, motif_targets)
```

## Arguments

- seurat_obj:

  A Seurat object

- motif_targets:

  list of motifs and their target genes

## Examples

``` r
SetMotifTargets
#> function (seurat_obj, motif_targets) 
#> {
#>     if (is.null(seurat_obj@misc$motifs)) {
#>         seurat_obj@misc$motifs <- list()
#>     }
#>     seurat_obj@misc$motifs$motif_targets <- motif_targets
#>     seurat_obj
#> }
#> <bytecode: 0x55d89745b9e0>
#> <environment: namespace:hdWGCNA>
```
