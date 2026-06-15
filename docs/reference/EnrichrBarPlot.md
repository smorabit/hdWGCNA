# EnrichrBarPlot

Generates bar plots from Enrichr data to visualize enriched terms for
hdWGCNA modules in a Seurat object. The function outputs a PDF file for
each module, with separate bar plots for each database.

## Usage

``` r
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots",
  n_terms = 25,
  p_cutoff = 0.05,
  p_adj = TRUE,
  plot_size = c(6, 15),
  logscale = FALSE,
  plot_bar_color = NULL,
  plot_text_color = NULL,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- outdir:

  A string specifying the directory where the output PDF files will be
  saved. Default is "enrichr_plots".

- n_terms:

  An integer indicating the number of top enriched terms to include in
  each bar plot. Default is 25.

- p_cutoff:

  A numeric value specifying the significance threshold for including
  terms (p-value or adjusted p-value). Only terms with p-values below
  this threshold are plotted. Default is 0.05.

- p_adj:

  Logical indicating whether to use the adjusted p-value (default: TRUE)
  or raw p-value for filtering terms.

- plot_size:

  A numeric vector of length 2 specifying the width and height of the
  output PDF files in inches. Default is c(6, 15).

- logscale:

  Logical specifying whether to log-transform the enrichment scores
  before plotting. Default is FALSE.

- plot_bar_color:

  A string specifying the color of the bars in the bar plots. If NULL
  (default), bars are colored according to the module's assigned color.

- plot_text_color:

  A string specifying the color of the text labels on the bar plots. If
  NULL (default), the color is automatically determined based on the bar
  color.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

Generates PDF files in the specified output directory, with one file per
WGCNA module. Each file contains bar plots for all databases associated
with that module.

## Details

This function processes the Enrichr output stored in a Seurat object,
filters enriched terms by significance, and generates bar plots for each
WGCNA module. Separate plots are created for each database included in
the Enrichr results. The top enriched terms for each module and database
are ordered by their combined enrichment score. If there are ties in the
enrichment scores, the function uses all tied terms.

The bar plot text and bar colors can be customized, and plots can be
saved with log-transformed enrichment scores if specified. Text wrapping
is applied to long term names for better readability.
