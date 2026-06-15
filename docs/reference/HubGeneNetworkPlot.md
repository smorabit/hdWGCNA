# HubGeneNetworkPlot

Construct a unified network plot comprising hub genes for multiple
modules.

## Usage

``` r
HubGeneNetworkPlot(
  seurat_obj,
  mods = "all",
  n_hubs = 6,
  n_other = 3,
  sample_edges = TRUE,
  edge_prop = 0.5,
  return_graph = FALSE,
  edge.alpha = 0.25,
  vertex.label.cex = 0.5,
  hub.vertex.size = 4,
  other.vertex.size = 1,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- mods:

  Names of the modules to plot. If mods = "all", all modules are
  plotted.

- n_hubs:

  The number of hub genes to plot for each module.

- n_other:

  The number of non-hub genes to sample from each module

- edge_prop:

  The proportion of edges in the graph to sample.

- return_graph:

  logical determining whether we return the graph (TRUE) or plot the
  graph (FALSE)

- edge.alpha:

  Scaling factor for the edge opacity

- vertex.label.cex:

  The font size of the gene labels

- hub.vertex.size:

  The size of the hub gene nodes

- other.vertex.size:

  The size of the other gene nodes

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Examples

``` r
HubGeneNetworkPlot
#> function (seurat_obj, mods = "all", n_hubs = 6, n_other = 3, 
#>     sample_edges = TRUE, edge_prop = 0.5, return_graph = FALSE, 
#>     edge.alpha = 0.25, vertex.label.cex = 0.5, hub.vertex.size = 4, 
#>     other.vertex.size = 1, wgcna_name = NULL, ...) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     MEs <- GetMEs(seurat_obj, wgcna_name)
#>     modules <- GetModules(seurat_obj, wgcna_name)
#>     if (all("all" %in% mods)) {
#>         mods <- levels(modules$module)
#>         mods <- mods[mods != "grey"]
#>     }
#>     else {
#>         if (!all(mods %in% unique(as.character(modules$module)))) {
#>             stop(paste0("Some selected modules are not found in wgcna_name: ", 
#>                 wgcna_name))
#>         }
#>         modules <- modules %>% subset(module %in% mods)
#>     }
#>     TOM <- GetTOM(seurat_obj, wgcna_name)
#>     hub_list <- lapply(mods, function(cur_mod) {
#>         cur <- subset(modules, module == cur_mod)
#>         cur[, c("gene_name", paste0("kME_", cur_mod))] %>% top_n(n_hubs) %>% 
#>             .$gene_name
#>     })
#>     names(hub_list) <- mods
#>     other_genes <- modules %>% subset(!(gene_name %in% unlist(hub_list))) %>% 
#>         group_by(module) %>% sample_n(n_other, replace = TRUE) %>% 
#>         .$gene_name %>% unique
#>     selected_genes <- c(unlist(hub_list), other_genes)
#>     selected_modules <- modules %>% subset(gene_name %in% selected_genes)
#>     subset_TOM <- TOM[selected_genes, selected_genes]
#>     selected_modules$geneset <- ifelse(selected_modules$gene_name %in% 
#>         other_genes, "other", "hub")
#>     selected_modules$size <- ifelse(selected_modules$geneset == 
#>         "hub", hub.vertex.size, other.vertex.size)
#>     selected_modules$label <- ifelse(selected_modules$geneset == 
#>         "hub", as.character(selected_modules$gene_name), "")
#>     selected_modules$fontcolor <- ifelse(selected_modules$color == 
#>         "black", "gray50", "black")
#>     edge_cutoff <- min(sapply(1:nrow(subset_TOM), function(i) {
#>         max(subset_TOM[i, ])
#>     }))
#>     edge_df <- reshape2::melt(subset_TOM) %>% subset(value >= 
#>         edge_cutoff)
#>     edge_df$color <- future.apply::future_sapply(1:nrow(edge_df), 
#>         function(i) {
#>             gene1 = as.character(edge_df[i, "Var1"])
#>             gene2 = as.character(edge_df[i, "Var2"])
#>             col1 <- modules %>% subset(gene_name == gene1) %>% 
#>                 .$color
#>             col2 <- modules %>% subset(gene_name == gene2) %>% 
#>                 .$color
#>             if (col1 == col2) {
#>                 col = col1
#>             }
#>             else {
#>                 col = "grey90"
#>             }
#>             col
#>         })
#>     groups <- unique(edge_df$color)
#>     print(groups)
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
#>     edge_df <- edge_df %>% group_by(color) %>% mutate(value = scale01(value))
#>     edge_df$color <- sapply(1:nrow(edge_df), function(i) {
#>         a = edge_df$value[i]
#>         alpha(edge_df$color[i], alpha = a)
#>     })
#>     g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, 
#>         vertices = selected_modules)
#>     l <- igraph::layout_with_fr(g, ...)
#>     if (return_graph) {
#>         return(g)
#>     }
#>     plot(g, layout = l, edge.color = adjustcolor(igraph::E(g)$color, 
#>         alpha.f = edge.alpha), vertex.size = igraph::V(g)$size, 
#>         edge.curved = 0, edge.width = 0.5, vertex.color = igraph::V(g)$color, 
#>         vertex.frame.color = igraph::V(g)$color, vertex.label = igraph::V(g)$label, 
#>         vertex.label.family = "Helvetica", vertex.label.font = 3, 
#>         vertex.label.color = igraph::V(g)$fontcolor, vertex.label.cex = vertex.label.cex, 
#>         ...)
#> }
#> <bytecode: 0x55d8a1f4ef60>
#> <environment: namespace:hdWGCNA>
```
