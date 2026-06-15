# ModuleCorrNetworks

Plot Module Eigengene correlogram

## Usage

``` r
ModuleCorrNetwork(
  seurat_obj,
  cluster_col = NULL,
  exclude_grey = TRUE,
  features = "hMEs",
  reduction = "umap",
  cor_cutoff = 0.2,
  label_vertices = FALSE,
  edge_scale = 5,
  vertex_size = 15,
  niter = 100,
  vertex_frame = FALSE,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

## Examples

``` r
ModuleCorrNetwork
#> function (seurat_obj, cluster_col = NULL, exclude_grey = TRUE, 
#>     features = "hMEs", reduction = "umap", cor_cutoff = 0.2, 
#>     label_vertices = FALSE, edge_scale = 5, vertex_size = 15, 
#>     niter = 100, vertex_frame = FALSE, wgcna_name = NULL) 
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
#>     modules <- GetModules(seurat_obj, wgcna_name)
#>     if (exclude_grey) {
#>         MEs <- MEs[, colnames(MEs) != "grey"]
#>         modules <- modules %>% subset(color != "grey")
#>     }
#>     mods <- colnames(MEs)
#>     if (is.null(cluster_col)) {
#>         clusters <- Idents(seurat_obj)
#>     }
#>     else {
#>         clusters <- droplevels(seurat_obj@meta.data[[cluster_col]])
#>     }
#>     MEs$cluster <- clusters
#>     cluster_ME_av <- do.call(rbind, lapply(split(MEs, MEs$cluster), 
#>         function(x) {
#>             colMeans(x[, mods])
#>         })) %>% as.data.frame
#>     top_clusters <- sapply(mods, function(x) {
#>         rownames(cluster_ME_av[cluster_ME_av[, x] == max(cluster_ME_av[, 
#>             x]), ])
#>     })
#>     red_df <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
#>     red_df$cluster <- clusters
#>     red_av <- do.call(rbind, lapply(split(red_df, red_df$cluster), 
#>         function(x) {
#>             colMeans(x[, 1:2])
#>         })) %>% as.data.frame
#>     cor_mat <- Hmisc::rcorr(as.matrix(MEs[, 1:ncol(MEs) - 1]))$r
#>     cor_mat[lower.tri(cor_mat)] <- NA
#>     cor_df <- reshape2::melt(cor_mat) %>% na.omit
#>     cor_df <- cor_df %>% subset(!(Var1 == Var2))
#>     cor_df <- cor_df %>% subset(abs(value) >= cor_cutoff)
#>     v_df <- data.frame(name = mods, cluster = as.character(top_clusters))
#>     unique_mods <- dplyr::distinct(modules[, c("module", "color")])
#>     rownames(unique_mods) <- unique_mods$module
#>     v_df$color <- unique_mods[v_df$name, "color"]
#>     v_df$x <- red_av[v_df$cluster, 1]
#>     v_df$y <- red_av[v_df$cluster, 2]
#>     g <- igraph::graph_from_data_frame(cor_df, directed = FALSE, 
#>         vertices = v_df)
#>     e <- get.edgelist(g, name = FALSE)
#>     l <- qgraph::qgraph.layout.fruchtermanreingold(e, vcount = vcount(g), 
#>         weights = igraph::E(g)$value, repulse.rad = (vcount(g)), 
#>         niter = niter, )
#>     plot_df <- rbind(cor_df, data.frame(Var1 = c("x", "y"), Var2 = c("y", 
#>         "x"), value = c(-1, 1)))
#>     temp <- ggplot(plot_df, aes(x = value, y = value, color = value)) + 
#>         geom_point() + scale_color_gradient2(high = "darkorchid1", 
#>         mid = "white", low = "seagreen", midpoint = 0)
#>     temp <- ggplot_build(temp)
#>     igraph::E(g)$color <- temp$data[[1]]$colour[1:nrow(cor_df)]
#>     if (label_vertices) {
#>         labels <- igraph::V(g)$name
#>     }
#>     else {
#>         labels <- NA
#>     }
#>     if (vertex_frame) {
#>         frame_color <- "black"
#>     }
#>     else {
#>         frame_color <- igraph::V(g)$color
#>     }
#>     plot(g, layout = l, edge.color = igraph::E(g)$color, edge.curved = 0, 
#>         edge.width = abs(igraph::E(g)$value) * edge_scale, vertex.color = igraph::V(g)$color, 
#>         vertex.frame.color = frame_color, vertex.label = labels, 
#>         vertex.label.family = "Helvetica", vertex.label.color = "black", 
#>         vertex.label.cex = 0.5, vertex.size = vertex_size)
#> }
#> <bytecode: 0x55d89a4a1510>
#> <environment: namespace:hdWGCNA>
```
