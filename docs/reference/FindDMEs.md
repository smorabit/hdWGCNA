# FindDMEs

Function to compare expression levels of co-expression modules between
two sets of cell barcodes.

## Usage

``` r
FindDMEs(
  seurat_obj,
  barcodes1,
  barcodes2,
  features = "MEs",
  harmonized = TRUE,
  wgcna_name = NULL,
  add_missing = FALSE,
  test.use = "wilcox",
  only.pos = FALSE,
  logfc.threshold = 0,
  min.pct = 0,
  verbose = FALSE,
  pseudocount.use = 0,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- barcodes1:

  character vector containing cell barcodes for the first group to test.
  Positive fold-change means up-regulated in this group.

- barcodes2:

  character vector containing cell barcodes for the second group to
  test. Negative fold-change means up-regulated in this group.

- features:

  indicate whether to use "MEs" or "ModuleScores" for the comparison

- harmonized:

  logical determining whether or not to use harmonized MEs

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- add_missing:

  logical determining whether or not to add missing modules back into
  the resulting dataframe with NA values.

- ...:

  Additional parameters for the FindMarkers function

## Value

A dataframe contaning differential ME results

## Examples

``` r
FindDMEs
#> function (seurat_obj, barcodes1, barcodes2, features = "MEs", 
#>     harmonized = TRUE, wgcna_name = NULL, add_missing = FALSE, 
#>     test.use = "wilcox", only.pos = FALSE, logfc.threshold = 0, 
#>     min.pct = 0, verbose = FALSE, pseudocount.use = 0, ...) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     CheckWGCNAName(seurat_obj, wgcna_name)
#>     if (!(all(barcodes1 %in% colnames(seurat_obj)))) {
#>         stop("Some barcodes in barcodes1 not found in colnames(seurat_obj).")
#>     }
#>     if (!(all(barcodes2 %in% colnames(seurat_obj)))) {
#>         stop("Some barcodes in barcodes2 not found in colnames(seurat_obj).")
#>     }
#>     if (length(intersect(barcodes1, barcodes2)) > 0) {
#>         stop("Some barcodes overlap in barcodes1 and barcodes2")
#>     }
#>     if (features == "MEs") {
#>         MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
#>     }
#>     else if (features == "ModuleScores") {
#>         MEs <- GetModuleScores(seurat_obj, wgcna_name)
#>     }
#>     else {
#>         stop("Invalid selection for features. Valid choices are MEs or ModuleScores.")
#>     }
#>     MEs <- MEs[, colnames(MEs) != "grey"]
#>     MEs[MEs < 0] <- 0
#>     MEs <- t(MEs)
#>     ME_assay <- Seurat::CreateAssayObject(MEs)
#>     DMEs <- FindMarkers(ME_assay, cells.1 = barcodes1, cells.2 = barcodes2, 
#>         slot = "counts", test.use = test.use, only.pos = only.pos, 
#>         logfc.threshold = logfc.threshold, min.pct = min.pct, 
#>         verbose = verbose, pseudocount.use = pseudocount.use, 
#>         ...)
#>     DMEs$module <- rownames(DMEs)
#>     if (add_missing) {
#>         missing_mods <- rownames(MEs)[!(rownames(MEs) %in% DMEs$module)]
#>         for (cur_mod in missing_mods) {
#>             DMEs[cur_mod, ] <- NA
#>             DMEs[cur_mod, "module"] <- cur_mod
#>         }
#>     }
#>     DMEs
#> }
#> <bytecode: 0x55d89df763a0>
#> <environment: namespace:hdWGCNA>
```
