# Displays the top n TFs in a set of modules as a bar plot

Displays the top n TFs in a set of modules as a bar plot

## Usage

``` r
MotifOverlapBarPlot(
  seurat_obj,
  n_tfs = 10,
  plot_size = c(5, 6),
  outdir = "MotifOverlaps/",
  motif_font = "helvetica_regular",
  module_names = NULL,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Examples

``` r
MotifOverlapBarPlot
#> function (seurat_obj, n_tfs = 10, plot_size = c(5, 6), outdir = "MotifOverlaps/", 
#>     motif_font = "helvetica_regular", module_names = NULL, wgcna_name = NULL) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     if (!dir.exists(outdir)) {
#>         dir.create(outdir)
#>     }
#>     modules <- GetModules(seurat_obj)
#>     mods <- levels(modules$module)
#>     mods <- mods[mods != "grey"]
#>     if (is.null(module_names)) {
#>         module_names <- mods
#>     }
#>     overlap_df <- GetMotifOverlap(seurat_obj, wgcna_name)
#>     motif_df <- GetMotifs(seurat_obj)
#>     pfm <- GetPFMList(seurat_obj)
#>     overlap_df$motif_ID <- motif_df$motif_ID[match(overlap_df$tf, 
#>         motif_df$motif_name)]
#>     overlap_df <- overlap_df %>% subset(module %in% module_names)
#>     for (cur_mod in module_names) {
#>         print(cur_mod)
#>         plot_df <- overlap_df %>% subset(module == cur_mod) %>% 
#>             top_n(n_tfs, wt = odds_ratio) %>% arrange(desc(odds_ratio))
#>         p1 <- plot_df %>% ggplot(aes(y = reorder(tf, odds_ratio), 
#>             fill = odds_ratio, x = odds_ratio)) + geom_bar(stat = "identity", 
#>             width = 0.7) + NoLegend() + scale_fill_gradient(high = unique(plot_df$color), 
#>             low = "grey90") + ylab("") + theme(axis.line.y = element_blank(), 
#>             axis.text.y = element_blank(), plot.margin = margin(t = 0, 
#>                 r = 0, b = 0, l = 0))
#>         plot_list <- list()
#>         for (i in 1:nrow(plot_df)) {
#>             cur_id <- plot_df[i, "motif_ID"]
#>             cur_name <- plot_df[i, "tf"]
#>             plot_list[[cur_id]] <- ggplot() + ggseqlogo::geom_logo(as.matrix(pfm[[cur_id]]), 
#>                 font = motif_font) + ggseqlogo::theme_logo() + 
#>                 xlab("") + ylab(cur_name) + theme(axis.text.x = element_blank(), 
#>                 axis.text.y = element_blank(), axis.title.y = element_text(angle = 0), 
#>                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
#>         }
#>         patch1 <- wrap_plots(plot_list, ncol = 1)
#>         outplot <- (patch1 | p1) + plot_layout(ncol = 2, widths = c(1, 
#>             2)) + plot_annotation(title = paste0("Motif overlaps with ", 
#>             cur_mod), theme = theme(plot.title = element_text(hjust = 0.5)))
#>         pdf(paste0(outdir, "/", cur_mod, "_motif_overlaps.pdf"), 
#>             width = plot_size[1], height = plot_size[2], useDingbats = FALSE)
#>         print(outplot)
#>         dev.off()
#>     }
#> }
#> <bytecode: 0x55d8a2e18570>
#> <environment: namespace:hdWGCNA>
```
