# EnrichrDotPlot

Generate a dot plot visualizing enrichment results from Enrichr for
hdWGCNA modules.

## Usage

``` r
EnrichrDotPlot(
  seurat_obj,
  database,
  mods = "all",
  n_terms = 3,
  p_cutoff = 0.05,
  p_adj = TRUE,
  break_ties = TRUE,
  term_size = 10,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- database:

  A character string specifying the name of the Enrichr database to use
  (e.g., "GO_Biological_Process_2021").

- mods:

  A character vector specifying the names of modules to include in the
  plot. If `mods = "all"` (default), all modules except the "grey"
  module are included.

- n_terms:

  An integer specifying the number of top enriched terms to plot for
  each module (default = 3).

- p_cutoff:

  A numeric value specifying the p-value threshold for filtering
  enriched terms (default = 0.05).

- p_adj:

  A logical value indicating whether to use adjusted p-values (`TRUE`,
  default) or raw p-values (`FALSE`) for filtering.

- break_ties:

  A logical value indicating whether to randomly select among tied terms
  to enforce `n_terms` (default = `TRUE`).

- term_size:

  A numeric value specifying the font size of the enriched terms
  displayed on the y-axis (default = 10).

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

A ggplot2 object representing the dot plot of enriched terms for the
specified modules and database.

## Details

This function creates dot plots from Enrichr results associated with
hdWGCNA modules. Each module is represented by its most enriched terms
from the specified Enrichr database. The size of the dots indicates the
enrichment score, and the color indicates the statistical significance
(-log10 transformed p-value).

- The function first retrieves WGCNA module and Enrichr data from the
  specified Seurat object.

- Modules are filtered based on the `mods` parameter, and enriched terms
  are filtered by significance using `p_cutoff` and `p_adj`.

- The top `n_terms` terms for each module are selected based on the
  Combined Score. If ties occur, the `break_ties` parameter determines
  how they are resolved.

- A dot plot is generated where each dot represents an enriched term,
  its size corresponds to the Combined Score (log-transformed), and its
  color indicates the significance (-log10 transformed p-value).
