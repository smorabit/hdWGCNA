# Module-Trait Correlation

Computes the correlation between WGCNA modules (eigengenes or scores)
and biological/clinical traits (metadata columns) for single cells.

## Usage

``` r
ModuleTraitCorrelation(
  seurat_obj,
  traits,
  group.by = NULL,
  features = "hMEs",
  cor_method = "pearson",
  subset_by = NULL,
  subset_groups = NULL,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object containing the hdWGCNA experiment.

- traits:

  A character vector of column names in `seurat_obj@meta.data` to
  correlate with each module. Traits must be numeric, integer, or
  factor. Character vectors should be converted to factors before
  running this function.

- group.by:

  A string containing the name of a column in `seurat_obj@meta.data`. If
  provided, correlations are computed separately for each group (e.g.,
  per cell type). If NULL (default), all cells are correlated together
  (global).

- features:

  Character; which feature to use to summarize each module?

  - "hMEs": Harmonized Module Eigengenes

  - "MEs": Standard Module Eigengenes.

  - "scores": Module Scores (Seurat `AddModuleScore`).

- cor_method:

  Character; method used for correlation. Valid choices: "pearson",
  "spearman", "kendall".

- subset_by:

  A string containing the name of a column to subset the data before
  correlation.

- subset_groups:

  A character vector specifying which groups from `subset_by` to include
  in the analysis.

- wgcna_name:

  The name of the hdWGCNA experiment in `seurat_obj@misc`. Default =
  NULL (uses active WGCNA).

- ...:

  Additional arguments passed to internal functions.

## Value

A Seurat object with the module-trait correlation results added to
`seurat_obj@misc[[wgcna_name]]$mt_cor`. The results include correlation
matrices, p-values, and FDR-corrected p-values for all cells and for
each group specified in `group.by`.

## Details

This function calculates the correlation between module eigengenes and
user-specified traits. If a trait is a factor, it is converted to
numeric based on the order of its levels. The function computes p-values
and False Discovery Rates (FDR) for each correlation.
