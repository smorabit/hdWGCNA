# Plots gene expression of hub genes as a heatmap

This function makes an expression heatmap of the top n hub genes per
module using Seurat's DoHeatmap, and then assembles them all into one
big heatmap.

## Usage

``` r
DoHubGeneHeatmap(
  seurat_obj,
  n_hubs = 10,
  n_cells = 200,
  group.by = NULL,
  module_names = NULL,
  combine = TRUE,
  draw.lines = TRUE,
  disp.min = -2.5,
  disp.max = 2.5,
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
DoHubGeneHeatmap
#> function (seurat_obj, n_hubs = 10, n_cells = 200, group.by = NULL, 
#>     module_names = NULL, combine = TRUE, draw.lines = TRUE, disp.min = -2.5, 
#>     disp.max = 2.5, wgcna_name = NULL) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     if (is.null(group.by)) {
#>         group.by <- "temp_ident"
#>         seurat_obj$temp_ident <- Idents(seurat_obj)
#>     }
#>     seurat_obj@meta.data[[group.by]] <- droplevels(seurat_obj@meta.data[[group.by]])
#>     modules <- GetModules(seurat_obj, wgcna_name)
#>     modules <- modules %>% subset(module != "grey") %>% mutate(module = droplevels(module))
#>     mods <- levels(modules$module)
#>     if (!is.null(module_names)) {
#>         print("here")
#>         mods <- module_names
#>         modules <- modules %>% subset(module %in% mods)
#>     }
#>     mod_colors <- modules %>% dplyr::select(c(module, color)) %>% 
#>         dplyr::distinct()
#>     hub_list <- lapply(mods, function(cur_mod) {
#>         cur <- subset(modules, module == cur_mod)
#>         cur <- cur[, c("gene_name", paste0("kME_", cur_mod))] %>% 
#>             top_n(n_hubs)
#>         colnames(cur)[2] <- "var"
#>         cur %>% arrange(desc(var)) %>% .$gene_name
#>     })
#>     names(hub_list) <- mods
#>     print(hub_list)
#>     seurat_obj$barcode <- colnames(seurat_obj)
#>     temp <- table(seurat_obj@meta.data[[group.by]])
#>     df <- data.frame()
#>     for (i in 1:length(temp)) {
#>         if (temp[[i]] < n_cells) {
#>             cur_df <- seurat_obj@meta.data %>% subset(get(group.by) == 
#>                 names(temp)[i])
#>         }
#>         else {
#>             cur_df <- seurat_obj@meta.data %>% subset(get(group.by) == 
#>                 names(temp)[i]) %>% sample_n(n_cells)
#>         }
#>         df <- rbind(df, cur_df)
#>     }
#>     seurat_plot <- seurat_obj %>% subset(barcode %in% df$barcode)
#>     plot_list <- list()
#>     for (i in 1:length(hub_list)) {
#>         print(i)
#>         cur_mod <- names(hub_list)[i]
#>         print(i)
#>         print(hub_list[[i]])
#>         print(i)
#>         if (i == 1) {
#>             plot_list[[i]] <- DoHeatmap(seurat_plot, features = hub_list[[i]], 
#>                 group.by = group.by, raster = TRUE, slot = "scale.data", 
#>                 disp.min = disp.min, disp.max = disp.max, label = FALSE, 
#>                 group.bar = FALSE, draw.lines = draw.lines)
#>         }
#>         else {
#>             plot_list[[i]] <- DoHeatmap(seurat_plot, features = hub_list[[i]], 
#>                 group.by = group.by, raster = TRUE, slot = "scale.data", 
#>                 group.bar.height = 0, label = FALSE, group.bar = FALSE, 
#>                 draw.lines = draw.lines, disp.min = disp.min, 
#>                 disp.max = disp.max) + NoLegend()
#>         }
#>         print(i)
#>         plot_list[[i]] <- plot_list[[i]] + theme(plot.margin = margin(0, 
#>             0, 0, 0), axis.text.y = element_text(face = "italic")) + 
#>             scale_y_discrete(position = "right")
#>         print(i)
#>     }
#>     n_total_cells <- ncol(seurat_plot)
#>     width_cbar <- n_total_cells/50
#>     mod_colors$value <- n_hubs
#>     mod_colors$dummy <- "colorbar"
#>     cbar_list <- list()
#>     for (i in 1:nrow(mod_colors)) {
#>         cbar_list[[i]] <- mod_colors[i, ] %>% ggplot(aes(y = value, 
#>             x = dummy)) + geom_bar(position = "stack", stat = "identity", 
#>             fill = mod_colors[i, ]$color) + umap_theme() + theme(plot.margin = margin(0, 
#>             0, 0, 0))
#>     }
#>     p_cbar <- wrap_plots(cbar_list, ncol = 1)
#>     if (combine) {
#>         out <- wrap_plots(plot_list, ncol = 1) + plot_layout(guides = "collect")
#>         out <- (p_cbar | out) + plot_layout(widths = c(width_cbar, 
#>             n_total_cells))
#>     }
#>     else {
#>         out <- plot_list
#>     }
#>     out
#> }
#> <bytecode: 0x55d895839f40>
#> <environment: namespace:hdWGCNA>
```
