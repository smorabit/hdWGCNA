# ModuleFeaturePlot

Plot module eigengenes as a FeaturePlot

## Usage

``` r
ModuleFeaturePlot(
  seurat_obj,
  module_names = NULL,
  wgcna_name = NULL,
  reduction = "umap",
  features = "hMEs",
  order_points = TRUE,
  restrict_range = TRUE,
  point_size = 0.5,
  alpha = 1,
  label_legend = FALSE,
  ucell = FALSE,
  raster = FALSE,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object

- features:

  What to plot? Can select hMEs, MEs, scores, or average

- order:

  TRUE, FALSE, or "shuffle" are valid options

## Examples

``` r
ModuleFeaturePlot
#> function (seurat_obj, module_names = NULL, wgcna_name = NULL, 
#>     reduction = "umap", features = "hMEs", order_points = TRUE, 
#>     restrict_range = TRUE, point_size = 0.5, alpha = 1, label_legend = FALSE, 
#>     ucell = FALSE, raster = FALSE, raster_dpi = 500, raster_scale = 1, 
#>     plot_ratio = 1, title = TRUE) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     if (features == "hMEs") {
#>         MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
#>     }
#>     else if (features == "MEs") {
#>         MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
#>     }
#>     else if (features == "scores") {
#>         MEs <- GetModuleScores(seurat_obj, wgcna_name)
#>     }
#>     else if (features == "average") {
#>         MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
#>         restrict_range <- FALSE
#>     }
#>     else (stop("Invalid feature selection. Valid choices: hMEs, MEs, scores, average"))
#>     if (ucell) {
#>         restrict_range <- FALSE
#>     }
#>     modules <- GetModules(seurat_obj, wgcna_name)
#>     if (is.null(module_names)) {
#>         module_names <- levels(modules$module)
#>         module_names <- module_names[module_names != "grey"]
#>     }
#>     umap <- seurat_obj@reductions[[reduction]]@cell.embeddings[, 
#>         1:2]
#>     x_name <- colnames(umap)[1]
#>     y_name <- colnames(umap)[2]
#>     plot_df <- cbind(umap, MEs) %>% as.data.frame()
#>     plot_list <- list()
#>     for (cur_mod in module_names) {
#>         print(cur_mod)
#>         cur_color <- modules %>% subset(module == cur_mod) %>% 
#>             .$color %>% unique
#>         plot_range <- plot_df[, cur_mod] %>% range
#>         if (restrict_range) {
#>             if (abs(plot_range[1]) > abs(plot_range[2])) {
#>                 plot_range[1] <- -1 * plot_range[2]
#>             }
#>             else {
#>                 plot_range[2] <- -1 * plot_range[1]
#>             }
#>             plot_df[, cur_mod] <- ifelse(plot_df[, cur_mod] > 
#>                 plot_range[2], plot_range[2], plot_df[, cur_mod])
#>             plot_df[, cur_mod] <- ifelse(plot_df[, cur_mod] < 
#>                 plot_range[1], plot_range[1], plot_df[, cur_mod])
#>         }
#>         cur_plot_df <- plot_df[, c(colnames(umap), cur_mod)]
#>         colnames(cur_plot_df)[3] <- "val"
#>         if (order_points == TRUE) {
#>             cur_plot_df <- cur_plot_df %>% dplyr::arrange(val)
#>         }
#>         else if (order_points == "shuffle") {
#>             cur_plot_df <- cur_plot_df[sample(nrow(cur_plot_df)), 
#>                 ]
#>         }
#>         p <- cur_plot_df %>% ggplot(aes_string(x = x_name, y = y_name, 
#>             color = "val"))
#>         if (raster) {
#>             p <- p + ggrastr::rasterise(geom_point(size = point_size, 
#>                 alpha = alpha), dpi = raster_dpi, scale = raster_scale)
#>         }
#>         else {
#>             p <- p + geom_point(size = point_size, alpha = alpha)
#>         }
#>         p <- p + umap_theme() + labs(color = "")
#>         if (title) {
#>             p <- p + ggtitle(cur_mod)
#>         }
#>         if (is.numeric(plot_ratio)) {
#>             p <- p + coord_fixed(ratio = plot_ratio)
#>         }
#>         if (!ucell) {
#>             p <- p + scale_color_gradient2(low = "grey75", mid = "grey95", 
#>                 high = cur_color, breaks = plot_range, labels = c("-", 
#>                   "+"), guide = guide_colorbar(ticks = FALSE, 
#>                   barwidth = 0.5, barheight = 4))
#>         }
#>         else {
#>             p <- p + scale_color_gradient(low = "grey95", high = cur_color, 
#>                 breaks = plot_range, labels = c("0", "+"), guide = guide_colorbar(ticks = FALSE, 
#>                   barwidth = 0.5, barheight = 4))
#>         }
#>         plot_list[[cur_mod]] <- p
#>     }
#>     if (length(plot_list) == 1) {
#>         p <- plot_list[[1]]
#>     }
#>     else {
#>         p <- plot_list
#>     }
#>     p
#> }
#> <bytecode: 0x55d8952993b0>
#> <environment: namespace:hdWGCNA>
```
