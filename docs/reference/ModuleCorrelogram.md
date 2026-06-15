# ModuleCorrelogram

Plot Module Eigengene correlogram

## Usage

``` r
ModuleCorrelogram(
  seurat_obj,
  MEs2 = NULL,
  features = "hMEs",
  order = "original",
  method = "ellipse",
  exclude_grey = TRUE,
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  sig.level = c(1e-04, 0.001, 0.01, 0.05),
  pch.cex = 0.7,
  col = NULL,
  ncolors = 200,
  wgcna_name = NULL,
  wgcna_name2 = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- features:

  What to plot? Can select hMEs, MEs, scores, or average

## Examples

``` r
MECorrelogram
#> Error: object 'MECorrelogram' not found
```
