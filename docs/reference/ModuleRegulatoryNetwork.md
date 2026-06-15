# ModuleRegulatoryNetwork

Summarizes Transcription Factor regulatory networks across co-expression
modules.

## Usage

``` r
ModuleRegulatoryNetwork(seurat_obj, TFs_only = TRUE, wgcna_name = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object containing the single-cell data and WGCNA results.

- TFs_only:

  Logical; if TRUE (default), only transcription factor (TF) genes are
  included in the regulatory network. If FALSE, the network includes all
  genes.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

A data frame with the following columns:

- source:

  The source module in the regulatory interaction.

- target:

  The target module in the regulatory interaction.

- n_pos:

  The number of positive regulatory links from the source module to the
  target module.

- sum_pos:

  The sum of the gain values for the positive links.

- n_neg:

  The number of negative regulatory links from the source module to the
  target module.

- sum_neg:

  The sum of the gain values for the negative links.

- mean_pos:

  The average gain of the positive links.

- mean_neg:

  The average gain of the negative links.

- score_pos:

  The number of positive links normalized by the number of TFs in the
  source module.

- score_neg:

  The number of negative links normalized by the number of TFs in the
  source module.

## Details

This function summarizes transcription factor regulatory networks across
modules to infer module-module relationships. The number of directed TF
links across modules are counted and normalized. The strength of these
links from the XGBoost model are also tracked.
