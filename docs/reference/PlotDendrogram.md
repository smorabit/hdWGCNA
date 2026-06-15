# PlotDendrogram

Plot WGCNA dendrogram

## Usage

``` r
PlotDendrogram(
  seurat_obj,
  groupLabels = "Module colors",
  wgcna_name = NULL,
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "",
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

## Examples

``` r
PlotDendrogram
#> function (seurat_obj, groupLabels = "Module colors", wgcna_name = NULL, 
#>     dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
#>     main = "", ...) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     net <- GetNetworkData(seurat_obj, wgcna_name)
#>     modules <- GetModules(seurat_obj, wgcna_name)
#>     WGCNA::plotDendroAndColors(net$dendrograms[[1]], as.character(modules$color), 
#>         groupLabels = groupLabels, dendroLabels = dendroLabels, 
#>         hang = hang, addGuide = addGuide, guideHang = guideHang, 
#>         main = main, ...)
#> }
#> <bytecode: 0x55d89910d1f0>
#> <environment: namespace:hdWGCNA>
```
