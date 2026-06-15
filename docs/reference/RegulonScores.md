# RegulonScores

Calculate expression scores for TF Regulons

## Usage

``` r
RegulonScores(
  seurat_obj,
  target_type = "positive",
  cor_thresh = 0.05,
  exclude_grey_genes = TRUE,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- target_type:

  Which types of TF target genes to compute scores for? "positive",
  "negative", or "both"

- cor_thresh:

  threshold for TF-gene correlation for genes to be included in the
  regulon score

- exclude_grey_genes:

  option to exclude genes that are in the grey module from the regulon
  scores

- wgcna_name:

  name of the WGCNA experiment

## Value

seurat_obj with the TF Regulon Scores added

## Details

RegulonScores calculates expression signatures for each TF regulon using
the UCell algorithm This function can calculate separate scores for TF
regulons for target genes that are positively or negatively correlated
with the TF, representing putative acivated or repressed genes. These
scores conveniently summarize the expression levels of the entire TF
regulon, similar to module eigengenes for the co-expression network
analyssis.
