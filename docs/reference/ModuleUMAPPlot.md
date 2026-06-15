# ModuleUMAPPlot

Makes a igraph network plot using the module UMAP

## Usage

``` r
ModuleUMAPPlot(
  seurat_obj,
  sample_edges = TRUE,
  edge_prop = 0.2,
  label_hubs = 5,
  edge.alpha = 0.25,
  vertex.label.cex = 0.5,
  label_genes = NULL,
  return_graph = FALSE,
  keep_grey_edges = TRUE,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- sample_edges:

  logical determining whether we downsample edges for plotting (TRUE),
  or take the strongst edges.

- edge_prop:

  proportion of edges to plot. If sample_edges=FALSE, the strongest
  edges are selected.

- label_hubs:

  the number of hub genes to label in each module

- edge.alpha:

  scaling factor for edge opacity

- vertex.label.cex:

  font size for labeled genes

- return_graph:

  logical determining whether to plot thr graph (FALSE) or return the
  igraph object (TRUE)

- keep_grey_edges:

  logical determining whether to show edges between genes in different
  modules (grey edges)

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Examples

``` r
ModuleUMAPPlot
#> function (seurat_obj, sample_edges = TRUE, edge_prop = 0.2, label_hubs = 5, 
#>     edge.alpha = 0.25, vertex.label.cex = 0.5, label_genes = NULL, 
#>     return_graph = FALSE, keep_grey_edges = TRUE, wgcna_name = NULL, 
#>     ...) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     TOM <- GetTOM(seurat_obj, wgcna_name)
#>     modules <- GetModules(seurat_obj, wgcna_name)
#>     umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
#>     mods <- levels(umap_df$module)
#>     mods <- mods[mods != "grey"]
#>     subset_TOM <- TOM[umap_df$gene, umap_df$gene[umap_df$hub == 
#>         "hub"]]
#>     hub_list <- lapply(mods, function(cur_mod) {
#>         cur <- subset(modules, module == cur_mod)
#>         cur[, c("gene_name", paste0("kME_", cur_mod))] %>% top_n(label_hubs) %>% 
#>             .$gene_name
#>     })
#>     names(hub_list) <- mods
#>     hub_labels <- as.character(unlist(hub_list))
#>     print("hub labels")
#>     print(hub_labels)
#>     print(label_genes)
#>     if (is.null(label_genes)) {
#>         label_genes <- hub_labels
#>     }
#>     else {
#>         if (!any(label_genes %in% umap_df$gene)) {
#>             stop("Some genes in label_genes not found in the UMAP.")
#>         }
#>         label_genes <- unique(c(label_genes, hub_labels))
#>     }
#>     print(label_genes)
#>     selected_modules <- modules[umap_df$gene, ]
#>     selected_modules <- cbind(selected_modules, umap_df[, c("UMAP1", 
#>         "UMAP2", "hub", "kME")])
#>     selected_modules$label <- ifelse(selected_modules$gene_name %in% 
#>         label_genes, selected_modules$gene_name, "")
#>     selected_modules$fontcolor <- ifelse(selected_modules$color == 
#>         "black", "gray50", "black")
#>     selected_modules$framecolor <- ifelse(selected_modules$gene_name %in% 
#>         label_genes, "black", selected_modules$color)
#>     edge_df <- subset_TOM %>% reshape2::melt()
#>     print(dim(edge_df))
#>     edge_df$color <- future.apply::future_sapply(1:nrow(edge_df), 
#>         function(i) {
#>             gene1 = as.character(edge_df[i, "Var1"])
#>             gene2 = as.character(edge_df[i, "Var2"])
#>             col1 <- selected_modules[selected_modules$gene_name == 
#>                 gene1, "color"]
#>             col2 <- selected_modules[selected_modules$gene_name == 
#>                 gene2, "color"]
#>             if (col1 == col2) {
#>                 col = col1
#>             }
#>             else {
#>                 col = "grey90"
#>             }
#>             col
#>         })
#>     if (!keep_grey_edges) {
#>         edge_df <- edge_df %>% subset(color != "grey90")
#>     }
#>     groups <- unique(edge_df$color)
#>     if (sample_edges) {
#>         temp <- do.call(rbind, lapply(groups, function(cur_group) {
#>             cur_df <- edge_df %>% subset(color == cur_group)
#>             n_edges <- nrow(cur_df)
#>             cur_sample <- sample(1:n_edges, round(n_edges * edge_prop))
#>             cur_df[cur_sample, ]
#>         }))
#>     }
#>     else {
#>         temp <- do.call(rbind, lapply(groups, function(cur_group) {
#>             cur_df <- edge_df %>% subset(color == cur_group)
#>             n_edges <- nrow(cur_df)
#>             cur_df %>% dplyr::top_n(round(n_edges * edge_prop), 
#>                 wt = value)
#>         }))
#>     }
#>     edge_df <- temp
#>     print(dim(edge_df))
#>     edge_df <- edge_df %>% group_by(color) %>% mutate(value = scale01(value))
#>     edge_df <- edge_df %>% arrange(value)
#>     edge_df <- rbind(subset(edge_df, color == "grey90"), subset(edge_df, 
#>         color != "grey90"))
#>     edge_df$color_alpha <- ifelse(edge_df$color == "grey90", 
#>         alpha(edge_df$color, alpha = edge_df$value/2), alpha(edge_df$color, 
#>             alpha = edge_df$value))
#>     selected_modules <- rbind(subset(selected_modules, hub == 
#>         "other"), subset(selected_modules, hub != "other"))
#>     selected_modules <- rbind(subset(selected_modules, label == 
#>         ""), subset(selected_modules, label != ""))
#>     g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, 
#>         vertices = selected_modules)
#>     if (return_graph) {
#>         return(g)
#>     }
#>     plot(g, layout = as.matrix(selected_modules[, c("UMAP1", 
#>         "UMAP2")]), edge.color = adjustcolor(igraph::E(g)$color_alpha, 
#>         alpha.f = edge.alpha), vertex.size = igraph::V(g)$kME * 
#>         3, edge.curved = 0, edge.width = 0.5, vertex.color = igraph::V(g)$color, 
#>         vertex.label = igraph::V(g)$label, vertex.label.dist = 1.1, 
#>         vertex.label.degree = -pi/4, vertex.label.family = "Helvetica", 
#>         vertex.label.font = 3, vertex.label.color = igraph::V(g)$fontcolor, 
#>         vertex.label.cex = 0, vertex.frame.color = igraph::V(g)$framecolor, 
#>         margin = 0)
#> }
#> <bytecode: 0x55d89fc3fdd8>
#> <environment: namespace:hdWGCNA>
```
