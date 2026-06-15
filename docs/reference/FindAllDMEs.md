# FindAllDMEs

Function to compare expression levels of co-expression modules between
two sets of cell barcodes.

## Usage

``` r
FindAllDMEs(
  seurat_obj,
  group.by,
  features = "MEs",
  harmonized = TRUE,
  add_missing = FALSE,
  wgcna_name = NULL,
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

- group.by:

  column in seurat_obj@meta.data containing cell grouping information

- features:

  indicate whether to use "MEs" or "ModuleScores" for the comparison

- harmonized:

  logical determining whether or not to use harmonized MEs

- add_missing:

  logical determining whether or not to add missing modules back into
  the resulting dataframe with NA values.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- ...:

  Additional parameters for the FindMarkers function

## Value

A dataframe contaning differential ME results

## Examples

``` r
FindAllDMEs
#> function (seurat_obj, group.by, features = "MEs", harmonized = TRUE, 
#>     add_missing = FALSE, wgcna_name = NULL, test.use = "wilcox", 
#>     only.pos = FALSE, logfc.threshold = 0, min.pct = 0, verbose = FALSE, 
#>     pseudocount.use = 0, ...) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     CheckWGCNAName(seurat_obj, wgcna_name)
#>     groups <- as.character(unique(seurat_obj@meta.data[[group.by]]))
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
#>     DMEs_list <- list()
#>     for (cur_group in groups) {
#>         print(cur_group)
#>         barcodes1 <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] == 
#>             cur_group]
#>         barcodes2 <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] != 
#>             cur_group]
#>         DMEs <- FindMarkers(ME_assay, cells.1 = barcodes1, cells.2 = barcodes2, 
#>             slot = "counts", test.use = test.use, only.pos = only.pos, 
#>             logfc.threshold = logfc.threshold, min.pct = min.pct, 
#>             verbose = verbose, pseudocount.use = pseudocount.use, 
#>             ...)
#>         DMEs$module <- rownames(DMEs)
#>         DMEs$group <- cur_group
#>         if (add_missing) {
#>             missing_mods <- rownames(MEs)[!(rownames(MEs) %in% 
#>                 DMEs$module)]
#>             for (cur_mod in missing_mods) {
#>                 DMEs[cur_mod, ] <- NA
#>                 DMEs[cur_mod, "module"] <- cur_mod
#>             }
#>         }
#>         DMEs_list[[cur_group]] <- DMEs
#>     }
#>     DMEs <- do.call(rbind, DMEs_list)
#>     DMEs
#> }
#> <bytecode: 0x55d89ad99098>
#> <environment: namespace:hdWGCNA>
```
