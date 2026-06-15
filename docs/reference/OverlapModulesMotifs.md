# Overlap modules with TF target genes

Overlap modules with TF target genes

## Usage

``` r
OverlapModulesMotifs(seurat_obj, wgcna_name = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Examples

``` r
OverlapModulesDEGs
#> function (seurat_obj, deg_df, fc_cutoff = 0.5, group_col = "cluster", 
#>     wgcna_name = NULL, ...) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     deg_df$group <- deg_df[, group_col]
#>     cell_groups <- deg_df$group %>% unique
#>     modules <- GetModules(seurat_obj)
#>     mods <- levels(modules$module)
#>     mods <- mods[mods != "grey"]
#>     if (fc_cutoff >= 0) {
#>         deg_df <- subset(deg_df, avg_log2FC >= fc_cutoff)
#>     }
#>     else {
#>         deg_df <- subset(deg_df, avg_log2FC <= fc_cutoff)
#>         deg_df$avg_log2FC <- -1 * deg_df$avg_log2FC
#>         fc_cutoff <- -1 * fc_cutoff
#>     }
#>     genome.size <- nrow(seurat_obj)
#>     overlap_df <- do.call(rbind, lapply(mods, function(cur_mod) {
#>         cur_module_genes <- modules %>% subset(module == cur_mod) %>% 
#>             .$gene_name
#>         cur_overlap_df <- do.call(rbind, lapply(cell_groups, 
#>             function(cur_group) {
#>                 cur_DEGs <- deg_df %>% subset(group == cur_group & 
#>                   p_val_adj <= 0.05 & avg_log2FC > fc_cutoff) %>% 
#>                   .$gene
#>                 cur_overlap <- testGeneOverlap(newGeneOverlap(cur_module_genes, 
#>                   cur_DEGs, genome.size = genome.size))
#>                 c(cur_overlap@odds.ratio, cur_overlap@pval, cur_overlap@Jaccard, 
#>                   length(cur_overlap@intersection))
#>             })) %>% as.data.frame
#>         colnames(cur_overlap_df) <- c("odds_ratio", "pval", "Jaccard", 
#>             "size_intersection")
#>         cur_overlap_df$module <- cur_mod
#>         cur_overlap_df$group <- cell_groups
#>         cur_overlap_df$color <- modules %>% subset(module == 
#>             cur_mod) %>% .$color %>% unique
#>         cur_overlap_df
#>     }))
#>     overlap_df$fdr <- p.adjust(overlap_df$pval, method = "fdr")
#>     overlap_df$Significance <- gtools::stars.pval(overlap_df$fdr)
#>     overlap_df$Significance <- ifelse(overlap_df$Significance == 
#>         ".", "", overlap_df$Significance)
#>     overlap_df$module <- factor(overlap_df$module, levels = mods)
#>     overlap_df <- overlap_df %>% dplyr::select(c(module, group, 
#>         color, odds_ratio, pval, fdr, Significance, Jaccard, 
#>         size_intersection))
#>     overlap_df
#> }
#> <bytecode: 0x55d8a3f7c9d8>
#> <environment: namespace:hdWGCNA>
```
