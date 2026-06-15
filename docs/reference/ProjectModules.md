# ProjectModules

Project a set of co-expression modules from a reference to a query
dataset

## Usage

``` r
ProjectModules(
  seurat_obj,
  seurat_ref = NULL,
  modules = NULL,
  group.by.vars = NULL,
  gene_mapping = NULL,
  genome1_col = NULL,
  genome2_col = NULL,
  overlap_proportion = 0.5,
  vars.to.regress = NULL,
  scale.model.use = "linear",
  wgcna_name = NULL,
  wgcna_name_proj = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object where the modules will be projected to

- seurat_ref:

  A Seurat object containing the co-expression modules to be projected

- modules:

  optionally provide a dataframe containing gene-module information to
  bypass seurat_ref

- group.by.vars:

  groups to harmonize by

- gene_mapping:

  dataframe to map gene names between genomes

- genome1_col:

  column in gene_mapping containing the genes for seurat_ref (reference)

- genome2_col:

  column in gene_mapping containing the genes for seurat_obj (query)

- overlap_proportion:

  the proportion of genes that must be present in seurat_obj for a given
  module to be projected. Default = 0.5 (50% of genes)

- vars.to.regress:

  character vector of variables in seurat_obj@meta.data to regress when
  running ScaleData

- scale.model.use:

  model to scale data when running ScaleData choices are "linear",
  "poisson", or "negbinom"

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_ref@misc slot

- wgcna_name_proj:

  The name of the hdWGCNA experiment to be created for the projected
  modules in seurat_obj

## Examples

``` r
ProjectModules
#> function (seurat_obj, seurat_ref = NULL, modules = NULL, group.by.vars = NULL, 
#>     gene_mapping = NULL, genome1_col = NULL, genome2_col = NULL, 
#>     overlap_proportion = 0.5, vars.to.regress = NULL, scale.model.use = "linear", 
#>     wgcna_name = NULL, wgcna_name_proj = NULL, ...) 
#> {
#>     if (!is.null(seurat_ref)) {
#>         if (is.null(wgcna_name)) {
#>             wgcna_name <- seurat_ref@misc$active_wgcna
#>         }
#>         CheckWGCNAName(seurat_ref, wgcna_name)
#>     }
#>     if (is.null(modules)) {
#>         modules <- GetModules(seurat_ref, wgcna_name)
#>     }
#>     else {
#>         if (!all(c("gene_name", "module", "color") %in% colnames(modules))) {
#>             stop("Missing columns in modules table. Required columns are gene_name, module, color")
#>         }
#>     }
#>     if (is.null(modules) & is.null(seurat_ref)) {
#>         stop("Either seurat_ref or modules must be provided.")
#>     }
#>     if (!is.factor(modules$module)) {
#>         modules$module <- as.factor(modules$module)
#>     }
#>     if (!is.null(gene_mapping)) {
#>         modules <- TransferModuleGenome(modules, gene_mapping, 
#>             genome1_col, genome2_col)
#>     }
#>     mods <- as.character(unique(modules$module))
#>     mod_props <- unlist(lapply(mods, function(x) {
#>         cur_mod <- subset(modules, module == x)
#>         sum(cur_mod$gene_name %in% rownames(seurat_obj))/nrow(cur_mod)
#>     }))
#>     mods_keep <- mods[mod_props >= overlap_proportion]
#>     mods_remove <- mods[mod_props < overlap_proportion]
#>     if (length(mods) > 0) {
#>         warning(paste0("The following modules will not be projected because too few genes are present in seurat_obj: ", 
#>             paste(mods_remove, collapse = ", ")))
#>     }
#>     modules <- subset(modules, module %in% mods_keep) %>% dplyr::mutate(module = droplevels(module))
#>     gene_names <- modules$gene_name
#>     genes_use <- intersect(gene_names, rownames(seurat_obj))
#>     modules <- modules %>% subset(gene_name %in% genes_use)
#>     if (!(wgcna_name_proj %in% names(seurat_obj@misc))) {
#>         seurat_obj <- SetupForWGCNA(seurat_obj, wgcna_name = wgcna_name_proj, 
#>             features = genes_use)
#>     }
#>     print("project modules")
#>     seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = group.by.vars, 
#>         modules = modules, vars.to.regress = vars.to.regress, 
#>         scale.model.use = scale.model.use, wgcna_name = wgcna_name_proj, 
#>         ...)
#>     seurat_obj
#> }
#> <bytecode: 0x55d89fb15e58>
#> <environment: namespace:hdWGCNA>
```
