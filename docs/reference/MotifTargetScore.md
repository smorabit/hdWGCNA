# MotifTargetScore

Computes gene expression scores for TF Motif target genes based on the
MotifScan.

## Usage

``` r
MotifTargetScore(
  seurat_obj,
  method = "Seurat",
  wgcna_genes = TRUE,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- method:

  Seurat or UCell?

## Examples

``` r
MotifTargetScore(pbmc)
#> Error in eval(expr, envir, enclos): object 'pbmc' not found
```
