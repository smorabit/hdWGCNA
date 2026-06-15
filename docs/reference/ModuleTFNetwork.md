# ModuleTFNetwork

Plotting the relationships between a TF and the co-expression modules

## Usage

``` r
ModuleTFNetwork(
  seurat_obj,
  tf_name,
  tf_gene_name,
  edge.alpha = 0.75,
  cor_thresh = 0.25,
  high_color = "red",
  mid_color = "grey",
  low_color = "blue",
  slot = "data",
  size.scale = 30,
  tf_x = 0,
  tf_y = 0,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- tf_name:

  the Motif name for the tf

- tf_gene_name:

  the gene associated with this tf in the rownames(seurat_obj)

- edge.alpha:

  scaling factor for edge opacity in the network

- cor_thresh:

  threshold to plot correlation edges between modules

- high_color:

  color for positive correlation

- mid_color:

  color for zero correlation

- low_color:

  color for negative correlation

- slot:

  the slot in the seurat object to extract expression data for the
  tf_gene_name

- size.scale:

  scaling factor for the size of each node

- tf_x:

  x coordinate for the TF if the TF is not found in the UMAP

- tf_y:

  y coordinate for the TF if the TF is not foudn in the UMAP

- wgcna_name:

  the name of the WGCNA experiment in the seurat object

## Examples

``` r
ModuleTFNetwork
#> function (seurat_obj, tf_name, tf_gene_name, edge.alpha = 0.75, 
#>     cor_thresh = 0.25, high_color = "red", mid_color = "grey", 
#>     low_color = "blue", slot = "data", size.scale = 30, tf_x = 0, 
#>     tf_y = 0, wgcna_name = NULL) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 
#>         "grey") %>% mutate(module = droplevels(module))
#>     MEs <- GetMEs(seurat_obj, TRUE, wgcna_name) %>% as.matrix
#>     MEs <- MEs[, colnames(MEs) != "grey"]
#>     mod_sizes <- table(modules$module)
#>     module_cor <- Hmisc::rcorr(x = MEs, type = "pearson")$r
#>     module_cor[lower.tri(module_cor)] <- NA
#>     module_cor <- reshape2::melt(module_cor) %>% na.omit
#>     module_cor <- subset(module_cor, abs(value) >= cor_thresh & 
#>         Var1 != Var2)
#>     module_cor
#>     cur_exp <- GetAssayData(seurat_obj, slot = slot)[tf_gene_name, 
#>         ]
#>     exp_cor <- Hmisc::rcorr(x = MEs, y = cur_exp)$r[1:ncol(MEs), 
#>         "y"]
#>     exp_cor <- data.frame(mod = names(exp_cor), value = as.numeric(exp_cor))
#>     plot_lim <- abs(max(c(abs(range(exp_cor$value)), abs(range(module_cor$value)))))
#>     p <- ggplot(module_cor, aes(x = Var1, y = Var2, color = value)) + 
#>         geom_point() + scale_color_gradient2(high = high_color, 
#>         mid = mid_color, low = low_color, limits = c(-1 * plot_lim, 
#>             plot_lim))
#>     ggp <- ggplot_build(p)
#>     module_cor$color <- ggp$data[[1]]$colour
#>     p <- ggplot(exp_cor, aes(x = mod, y = mod, color = value)) + 
#>         geom_point() + scale_color_gradient2(high = high_color, 
#>         mid = mid_color, low = low_color, limits = c(-1 * plot_lim, 
#>             plot_lim))
#>     ggp <- ggplot_build(p)
#>     exp_cor$color <- ggp$data[[1]]$colour
#>     umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
#>     mods <- levels(umap_df$modules)
#>     tf_match <- GetMotifMatrix(seurat_obj)
#>     tf_targets <- GetMotifTargets(seurat_obj)
#>     motif_df <- GetMotifs(seurat_obj)
#>     overlap_df <- GetMotifOverlap(seurat_obj)
#>     centroid_df <- umap_df %>% dplyr::select(c(UMAP1, UMAP2, 
#>         module)) %>% dplyr::group_by(module) %>% dplyr::summarise(x = mean(UMAP1), 
#>         y = mean(UMAP2))
#>     cur_overlap <- subset(overlap_df, tf == tf_name)
#>     node_df <- dplyr::left_join(centroid_df, cur_overlap, by = "module") %>% 
#>         dplyr::rename(c(UMAP1 = x, UMAP2 = y, name = module))
#>     node_df$size <- as.numeric(node_df$size_intersection)/as.numeric(mod_sizes)
#>     if (tf_gene_name %in% umap_df$gene) {
#>         tf_df <- umap_df[umap_df$gene == tf_gene_name, ] %>% 
#>             dplyr::select(-c(hub, kME, module)) %>% dplyr::rename(name = gene)
#>     }
#>     else {
#>         tf_df <- data.frame(name = tf_gene_name, module = "grey", 
#>             color = "grey", UMAP1 = tf_x, UMAP2 = tf_y)
#>     }
#>     tf_df$size <- 0.25
#>     edge_df <- data.frame(Var1 = tf_gene_name, Var2 = as.character(node_df$name), 
#>         value = node_df$odds_ratio, color = exp_cor$color)
#>     edge_df$color <- ifelse(node_df$fdr <= 0.05, edge_df$color, 
#>         "grey")
#>     edge_df <- subset(edge_df, color != "grey")
#>     node_df <- dplyr::bind_rows(node_df, tf_df) %>% as.data.frame()
#>     rownames(node_df) <- as.character(node_df$name)
#>     g1 <- igraph::graph_from_data_frame(edge_df, directed = TRUE, 
#>         vertices = node_df)
#>     g2 <- igraph::graph_from_data_frame(module_cor, directed = FALSE, 
#>         vertices = node_df)
#>     plot(g2, layout = as.matrix(node_df[, c("UMAP1", "UMAP2")]), 
#>         vertex.size = 1, edge.curved = 0, edge.width = 1, vertex.color = "grey", 
#>         vertex.label = "", edge.color = adjustcolor(igraph::E(g2)$color, 
#>             alpha.f = edge.alpha), )
#>     plot(g1, layout = as.matrix(node_df[, c("UMAP1", "UMAP2")]), 
#>         edge.color = adjustcolor(igraph::E(g1)$color), vertex.size = igraph::V(g1)$size * 
#>             size.scale, edge.curved = 0, edge.width = edge_df$value * 
#>             2, vertex.color = igraph::V(g1)$color, vertex.label = igraph::V(g1)$name, 
#>         vertex.label.dist = 1.1, vertex.label.degree = -pi/4, 
#>         vertex.label.family = "Helvetica", vertex.label.font = 3, 
#>         vertex.label.color = "black", vertex.label.cex = 0, vertex.frame.color = "black", 
#>         margin = 0, edge.arrow.size = edge_df$value/2, add = TRUE)
#> }
#> <bytecode: 0x7fced738ad68>
#> <environment: namespace:hdWGCNA>
```
