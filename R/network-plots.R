#' ModuleNetworkPlot
#'
#' Visualizes the top hub genes of each module as a circular network: the
#' strongest hubs on an inner circle, the next ones on an outer circle, and the
#' strongest TOM connections between them as edges.
#'
#' @param seurat_obj A Seurat object
#' @param n_inner Number of genes to plot on the inner circle
#' @param n_outer Number of genes to plot on the outer circle.
#' @param n_conns Number of gene-gene co-expression connections to plot, sorted by the strongest connections.
#' @param mods Names of the modules to plot. If mods = "all", all modules are plotted.
#' @param outdir Deprecated. If given, the plots are also written to this directory as one PDF per module (use \code{SavePlotList} instead). Default = NULL
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param plot_size Width and height of the PDFs written to \code{outdir}. Default = c(6, 6)
#' @param edge.alpha Value between 0 and 1 scaling the opacity of the edges.
#' @param edge.width Width of the edges.
#' @param vertex.label.cex Size of the gene labels, relative to the default.
#' @param vertex.size Size of the gene nodes.
#' @param ... Currently unused.
#' @return A named list of ggplot objects, one per module (modules with too few genes are skipped).
#' @export
ModuleNetworkPlot <- function(
  seurat_obj,
  n_inner = 10,
  n_outer = 15,
  n_conns = 500,
  mods = 'all',
  outdir = NULL,
  wgcna_name = NULL,
  plot_size = c(6, 6),
  edge.alpha = 0.25,
  edge.width = 1,
  vertex.label.cex = 1,
  vertex.size = 6, ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- select_modules(modules, mods, wgcna_name)

  # check if we have eigengene-based connectivities:
  if (!all(paste0('kME_', as.character(mods)) %in% colnames(modules))) {
    stop('Eigengene-based connectivity (kME) not found. Did you run ModuleEigengenes and ModuleConnectivity?')
  }

  TOM <- GetTOM(seurat_obj, wgcna_name)
  hub_df <- GetHubGenes(seurat_obj, n_hubs = n_inner + n_outer, wgcna_name = wgcna_name)

  plot_list <- list()
  for (cur_mod in mods) {
    cur_color <- unique(as.character(modules$color[modules$module == cur_mod]))
    cur_genes <- as.character(hub_df$gene_name[hub_df$module == cur_mod])
    n_genes <- length(cur_genes)

    if (n_genes < (n_inner + 1)) {
      message('Skipping ', cur_mod, ', too few genes to plot.')
      next
    }

    # keep the strongest connections among the hub genes
    reducedTOM <- as.matrix(TOM[cur_genes, cur_genes])
    diag(reducedTOM) <- 0
    reducedTOM[lower.tri(reducedTOM)] <- 0
    orderind <- order(reducedTOM, decreasing = TRUE)
    cur_n_conns <- min(n_conns, sum(reducedTOM > 0))
    reducedTOM[orderind[-seq_len(cur_n_conns)]] <- 0
    reducedTOM <- scale01(reducedTOM)

    edge_df <- reshape2::melt(reducedTOM, varnames = c('from', 'to'))
    edge_df <- edge_df[edge_df$value > 0, ]
    edge_df <- edge_df[order(edge_df$value), ]

    # hubs on an inner circle of half the radius, the rest on an outer circle
    n_out <- n_genes - n_inner
    node_df <- data.frame(
      name = cur_genes,
      x = c(cos(2 * pi * (seq_len(n_inner) - 1) / n_inner) / 2, cos(2 * pi * (seq_len(n_out) - 1) / n_out)),
      y = c(sin(2 * pi * (seq_len(n_inner) - 1) / n_inner) / 2, sin(2 * pi * (seq_len(n_out) - 1) / n_out))
    )

    g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = node_df)
    p <- ggraph::ggraph(g, layout = 'manual', x = node_df$x, y = node_df$y) +
      ggraph::geom_edge_link(
        ggplot2::aes(alpha = .data$value * edge.alpha),
        colour = cur_color, edge_width = edge.width * 0.5
      ) +
      ggraph::geom_node_point(fill = cur_color, shape = 21, color = 'black', size = vertex.size * 0.6) +
      ggraph::geom_node_text(
        ggplot2::aes(label = .data$name),
        repel = TRUE, fontface = 'italic', size = 3.5 * vertex.label.cex
      ) +
      ggraph::scale_edge_alpha_identity() +
      ggplot2::coord_fixed(clip = 'off') +
      ggplot2::labs(title = cur_mod) +
      ggraph::theme_graph(base_family = '') +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    plot_list[[cur_mod]] <- p
  }

  if (!is.null(outdir)) {
    deprecated_outdir('ModuleNetworkPlot')
    SavePlotList(plot_list, outdir = outdir, width = plot_size[1], height = plot_size[2])
  }
  plot_list
}


#' HubGeneNetworkPlot
#'
#' Construct a unified network plot comprising hub genes for multiple modules.
#'
#' @param seurat_obj A Seurat object
#' @param mods Names of the modules to plot. If mods = "all", all modules are plotted.
#' @param n_hubs The number of hub genes to plot for each module.
#' @param n_other The number of non-hub genes to sample from each module
#' @param sample_edges logical determining whether we randomly downsample edges (TRUE), or take the strongest edges (FALSE).
#' @param edge_prop The proportion of edges in the graph to sample.
#' @param return_graph logical determining whether we return the graph (TRUE) or the plot (FALSE)
#' @param edge.alpha Scaling factor for the edge opacity
#' @param vertex.label.cex Size of the hub gene labels, relative to the default.
#' @param hub.vertex.size The size of the hub gene nodes
#' @param other.vertex.size The size of the other gene nodes
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Additional arguments passed to \code{igraph::layout_with_fr}.
#' @return A ggplot object, or an igraph object if \code{return_graph = TRUE}.
#' @keywords scRNA-seq
#' @export
HubGeneNetworkPlot <- function(
  seurat_obj, mods = 'all',
  n_hubs = 6, n_other = 3,
  sample_edges = TRUE,
  edge_prop = 0.5,
  return_graph = FALSE,
  edge.alpha = 0.25,
  vertex.label.cex = 0.5,
  hub.vertex.size = 4,
  other.vertex.size = 1,
  wgcna_name = NULL,
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- select_modules(modules, mods, wgcna_name)
  modules <- modules[modules$module %in% mods, ]

  TOM <- GetTOM(seurat_obj, wgcna_name)

  # hub genes, plus a sample of other genes from each module
  hub_genes <- as.character(GetHubGenes(seurat_obj, n_hubs = n_hubs, mods = mods, wgcna_name = wgcna_name)$gene_name)
  other_genes <- unique(unlist(lapply(mods, function(cur_mod) {
    cur <- setdiff(modules$gene_name[modules$module == cur_mod], hub_genes)
    cur[sample.int(length(cur), min(n_other, length(cur)))]
  })))

  selected_genes <- c(hub_genes, other_genes)
  node_df <- modules[match(selected_genes, modules$gene_name), c('gene_name', 'module', 'color')]
  colnames(node_df)[1] <- 'name'
  node_df$color <- as.character(node_df$color)
  node_df$hub <- node_df$name %in% hub_genes
  node_df$size <- ifelse(node_df$hub, hub.vertex.size, other.vertex.size)
  node_df$label <- ifelse(node_df$hub, node_df$name, '')

  # keep edges at least as strong as the weakest "strongest edge" of any gene, so that
  # every gene keeps at least one edge
  subset_TOM <- as.matrix(TOM[selected_genes, selected_genes])
  diag(subset_TOM) <- 0
  edge_cutoff <- min(apply(subset_TOM, 1, max))
  subset_TOM[lower.tri(subset_TOM, diag = TRUE)] <- NA
  edge_df <- reshape2::melt(subset_TOM)
  edge_df <- edge_df[!is.na(edge_df$value) & edge_df$value >= edge_cutoff, ]
  edge_df$color <- edge_colors(edge_df, modules)

  edge_df <- subsample_edges(edge_df, edge_prop, sample_edges)
  colnames(edge_df)[1:2] <- c('from', 'to')

  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = node_df)
  if (return_graph) {
    return(g)
  }
  layout <- igraph::layout_with_fr(g, ...)

  ggraph::ggraph(g, layout = 'manual', x = layout[, 1], y = layout[, 2]) +
    ggraph::geom_edge_link(
      ggplot2::aes(colour = .data$color, alpha = .data$value * edge.alpha),
      edge_width = 0.3
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(color = .data$color, size = .data$size)
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$label),
      repel = TRUE, fontface = 'italic', size = 7 * vertex.label.cex, max.overlaps = Inf
    ) +
    ggraph::scale_edge_colour_identity() +
    ggraph::scale_edge_alpha_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_identity() +
    ggraph::theme_graph(base_family = '')
}


#' RunModuleUMAP
#'
#' Run UMAP on co-expression matrix using hub genes as features.
#'
#' @param seurat_obj A Seurat object
#' @param n_hubs number of hub genes to use in the UMAP computation
#' @param exclude_grey logical indicating whether to include grey module
#' @param genes_use character vector of genes to use for the UMAP, must already be present in GetModules(seurat_obj)
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param n_neighbors Number of nearest neighbors passed to \code{uwot::umap}. Default = 25.
#' @param metric Distance metric passed to \code{uwot::umap}. Default = 'cosine'.
#' @param spread Effective scale of embedded points passed to \code{uwot::umap}. Default = 1.
#' @param min_dist Minimum distance between embedded points passed to \code{uwot::umap}. Default = 0.4.
#' @param supervised logical determining whether to run a supervised UMAP using the module assignments as labels.
#' @param ... Additional parameters supplied to uwot::umap
#' @return The Seurat object with the module UMAP coordinates stored in \code{seurat_obj@misc[[wgcna_name]]}; retrieve them with \code{GetModuleUMAP}.
#' @export
RunModuleUMAP <- function(
  seurat_obj,
  n_hubs = 10,
  exclude_grey = TRUE,
  genes_use = NULL,
  wgcna_name = NULL,
  n_neighbors = 25,
  metric = 'cosine',
  spread = 1,
  min_dist = 0.4,
  supervised = FALSE,
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get modules,
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)

  # check if we have eigengene-based connectivities:
  if (!all(paste0('kME_', as.character(mods)) %in% colnames(modules))) {
    stop('Eigengene-based connectivity (kME) not found. Did you run ModuleEigengenes and ModuleConnectivity?')
  }

  # handle excluding grey module
  if (exclude_grey) {
    mods <- mods[mods != 'grey']
    modules <- subset(modules, module != 'grey')
  }

  # get hub genes:
  hub_df <- GetHubGenes(seurat_obj, n_hubs = n_hubs, wgcna_name = wgcna_name)

  # get all genes that aren't in gray mod
  selected_genes <- modules[modules$module %in% mods, 'gene_name']

  if (!is.null(genes_use)) {
    selected_genes <- selected_genes[selected_genes %in% genes_use]
  }

  # get the TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # subset the TOM for umap
  # keep all genes as rows, and keep only hubs as cols (features)
  feature_mat <- TOM[selected_genes, hub_df$gene_name]

  # run UMAP
  if (supervised) {
    hub_umap <- uwot::umap(
      X = feature_mat,
      min_dist = min_dist,
      n_neighbors = n_neighbors,
      metric = metric,
      spread = spread,
      y = modules$module, # for supervised UMAP
      ...
    )
  } else {
    hub_umap <- uwot::umap(
      X = feature_mat,
      min_dist = min_dist,
      n_neighbors = n_neighbors,
      metric = metric,
      spread = spread,
      ...
    )
  }

  # set up plotting df
  plot_df <- as.data.frame(hub_umap)
  colnames(plot_df) <- c('UMAP1', 'UMAP2')
  plot_df$gene <- rownames(feature_mat)

  # add module color, and hub gene status to the plotting df:
  ix <- match(plot_df$gene, modules$gene_name)
  plot_df$module <- modules$module[ix]
  plot_df$color <- modules$color[ix]
  plot_df$hub <- ifelse(
    plot_df$gene %in% as.character(hub_df$gene_name), 'hub', 'other'
  )

  # get kME values for each gene
  kMEs <- do.call(rbind, lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur <- cur[, c('gene_name', paste0('kME_', cur_mod))]
    colnames(cur) <- c('gene_name', 'kME')

    # scale kMEs between 0 & 1:
    cur$kME <- scale01(cur$kME)
    cur
  }))
  ix <- kMEs$gene_name[match(plot_df$gene, kMEs$gene_name)]
  plot_df$kME <- kMEs[ix, 'kME']

  # add the UMAP to the Seurat object:
  seurat_obj <- SetModuleUMAP(seurat_obj, plot_df, wgcna_name)

  seurat_obj
}


#' ModuleUMAPPlot
#'
#' Plots the module UMAP from \code{RunModuleUMAP} as a network: genes are
#' placed at their UMAP coordinates and connected to hub genes by their TOM.
#'
#' @param seurat_obj A Seurat object
#' @param sample_edges logical determining whether we downsample edges for plotting (TRUE), or take the strongest edges.
#' @param edge_prop proportion of edges to plot. If sample_edges=FALSE, the strongest edges are selected.
#' @param label_hubs the number of hub genes to label in each module
#' @param edge.alpha scaling factor for edge opacity
#' @param vertex.label.cex Size of the gene labels, relative to the default.
#' @param label_genes Optional character vector of additional genes to label, alongside the top hub genes.
#' @param return_graph logical determining whether to return the plot (FALSE) or the igraph object (TRUE)
#' @param keep_grey_edges logical determining whether to show edges between genes in different modules (grey edges)
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Currently unused.
#' @return A ggplot object, or an igraph object if \code{return_graph = TRUE}.
#' @keywords scRNA-seq
#' @export
ModuleUMAPPlot <- function(
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
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  TOM <- GetTOM(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  if (is.null(umap_df)) {
    stop('Module UMAP not found. Run RunModuleUMAP first.')
  }
  mods <- setdiff(levels(umap_df$module), 'grey')

  # genes to label: top hub genes of each module, plus any requested genes
  hub_labels <- as.character(GetHubGenes(seurat_obj, n_hubs = label_hubs, mods = mods, wgcna_name = wgcna_name)$gene_name)
  if (!is.null(label_genes)) {
    if (!all(label_genes %in% umap_df$gene)) {
      stop(paste0(
        'Genes in label_genes not found in the module UMAP: ',
        paste(setdiff(label_genes, umap_df$gene), collapse = ', ')
      ))
    }
  }
  label_genes <- unique(c(label_genes, hub_labels))

  node_df <- data.frame(
    name = umap_df$gene,
    x = umap_df$UMAP1,
    y = umap_df$UMAP2,
    color = as.character(umap_df$color),
    kME = umap_df$kME,
    hub = umap_df$hub
  )
  node_df$label <- ifelse(node_df$name %in% label_genes, node_df$name, '')
  node_df$frame <- ifelse(node_df$label != '', 'black', node_df$color)
  # draw hubs and labeled genes on top
  node_df <- node_df[order(node_df$hub == 'hub', node_df$label != ''), ]

  # edges between every gene and the hub genes
  hubs <- umap_df$gene[umap_df$hub == 'hub']
  edge_df <- reshape2::melt(as.matrix(TOM[umap_df$gene, hubs]))
  edge_df <- edge_df[as.character(edge_df$Var1) != as.character(edge_df$Var2), ]
  edge_df$color <- edge_colors(edge_df, modules)
  if (!keep_grey_edges) {
    edge_df <- edge_df[edge_df$color != 'grey90', ]
  }
  edge_df <- subsample_edges(edge_df, edge_prop, sample_edges)
  # grey edges are fainter and drawn underneath
  edge_df$alpha <- edge_df$value * edge.alpha * ifelse(edge_df$color == 'grey90', 0.5, 1)
  edge_df <- edge_df[order(edge_df$color != 'grey90', edge_df$value), ]
  colnames(edge_df)[1:2] <- c('from', 'to')

  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = node_df)
  if (return_graph) {
    return(g)
  }

  ggraph::ggraph(g, layout = 'manual', x = node_df$x, y = node_df$y) +
    ggraph::geom_edge_link(
      ggplot2::aes(colour = .data$color, alpha = .data$alpha),
      edge_width = 0.3
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(fill = .data$color, color = .data$frame, size = .data$kME),
      shape = 21, stroke = 0.3
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$label),
      repel = TRUE, fontface = 'italic', size = 7 * vertex.label.cex, max.overlaps = Inf
    ) +
    ggraph::scale_edge_colour_identity() +
    ggraph::scale_edge_alpha_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_continuous(range = c(0.3, 3), guide = 'none') +
    ggraph::theme_graph(base_family = '')
}


#' ModuleCorrNetwork
#'
#' Plot the correlations between modules as a network: modules are nodes and
#' edges connect modules whose features are correlated, colored by the
#' correlation.
#'
#' @param seurat_obj A Seurat object
#' @param cluster_col Column in \code{seurat_obj@meta.data} with cluster labels. Each module starts the layout at the centroid (in \code{reduction}) of the cluster where it is most highly expressed. Defaults to \code{Idents(seurat_obj)}.
#' @param exclude_grey logical determining whether to exclude the grey module.
#' @param features Which module features to correlate: 'hMEs', 'MEs', 'scores', or 'average'.
#' @param reduction Name of the dimensionality reduction used to compute cluster centroids. Default = 'umap'. Set to NULL to start the layout from random positions.
#' @param cor_cutoff Minimum absolute correlation for an edge to be drawn. Default = 0.2.
#' @param label_vertices logical determining whether to label the module nodes.
#' @param edge_scale Scaling factor applied to the absolute correlation to set edge width.
#' @param vertex_size Size of the module nodes.
#' @param niter Number of iterations for the Fruchterman-Reingold layout.
#' @param vertex_frame logical determining whether to draw a black outline around the nodes (TRUE) or use the module color (FALSE).
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A ggplot object.
#' @keywords scRNA-seq
#' @export
ModuleCorrNetwork <- function(
  seurat_obj, cluster_col = NULL, exclude_grey = TRUE,
  features = 'hMEs',
  reduction = 'umap', cor_cutoff = 0.2, label_vertices = FALSE, edge_scale = 5,
  vertex_size = 15, niter = 100, vertex_frame = FALSE, wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  MEs <- switch(features,
    hMEs = GetMEs(seurat_obj, TRUE, wgcna_name),
    MEs = GetMEs(seurat_obj, FALSE, wgcna_name),
    scores = GetModuleScores(seurat_obj, wgcna_name),
    average = GetAvgModuleExpr(seurat_obj, wgcna_name),
    stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
  )
  MEs <- as.matrix(MEs)
  if (exclude_grey) {
    MEs <- MEs[, colnames(MEs) != 'grey', drop = FALSE]
  }
  mods <- colnames(MEs)

  modules <- GetModules(seurat_obj, wgcna_name)
  mod_colors <- as.character(modules$color[match(mods, modules$module)])

  # correlation edges
  cor_mat <- stats::cor(MEs)
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  edge_df <- reshape2::melt(cor_mat, varnames = c('from', 'to'))
  edge_df <- edge_df[!is.na(edge_df$value) & abs(edge_df$value) >= cor_cutoff, ]

  node_df <- data.frame(name = mods, color = mod_colors)
  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = node_df)

  # start each module at the centroid of the cluster where it is highest
  coords <- NULL
  if (!is.null(reduction)) {
    clusters <- if (is.null(cluster_col)) Idents(seurat_obj) else seurat_obj@meta.data[[cluster_col]]
    clusters <- factor(clusters)
    check_reduction(seurat_obj, reduction)
    emb <- seurat_obj@reductions[[reduction]]@cell.embeddings[, 1:2]
    clusters <- droplevels(clusters)
    n_cells <- as.vector(table(clusters))
    cluster_ME <- rowsum(MEs, clusters) / n_cells
    top_cluster <- rownames(cluster_ME)[apply(cluster_ME, 2, which.max)]
    centroids <- rowsum(emb, clusters) / n_cells
    coords <- centroids[top_cluster, , drop = FALSE]
    # modules peaking in the same cluster would start on top of each other
    coords <- coords + matrix(stats::rnorm(length(coords), sd = 1e-3 * max(abs(coords))), ncol = 2)
  }
  layout <- igraph::layout_with_fr(
    g,
    coords = coords,
    weights = abs(igraph::E(g)$value),
    niter = niter
  )

  frame <- if (vertex_frame) 'black' else mod_colors
  p <- ggraph::ggraph(g, layout = 'manual', x = layout[, 1], y = layout[, 2]) +
    ggraph::geom_edge_link(ggplot2::aes(colour = .data$value, edge_width = abs(.data$value) * edge_scale)) +
    ggraph::geom_node_point(fill = mod_colors, color = frame, shape = 21, size = vertex_size / 3) +
    ggraph::scale_edge_colour_gradient2(
      low = 'seagreen', mid = 'white', high = 'darkorchid1', midpoint = 0, limits = c(-1, 1),
      name = 'Correlation', guide = ggraph::guide_edge_colourbar()
    ) +
    ggraph::scale_edge_width_identity() +
    ggraph::theme_graph(base_family = '')
  if (label_vertices) {
    p <- p + ggraph::geom_node_text(ggplot2::aes(label = .data$name), repel = TRUE, size = 3)
  }
  p
}


#' SavePlotList
#'
#' Save a list of plots, either as one file per plot in a directory or as a
#' single multi-page PDF.
#'
#' @param plot_list A list of ggplot objects, such as the output of
#'   \code{ModuleNetworkPlot} or \code{EnrichrBarPlot}. Names are used as file names.
#' @param outdir Directory to write one file per plot, named after the list
#'   names. Created if it does not exist.
#' @param file Path of a single multi-page PDF to write instead of \code{outdir}.
#' @param width,height Plot size in inches. Default = 6
#' @param device File type for \code{outdir}, passed to \code{ggplot2::ggsave}. Default = 'pdf'
#' @param ... Additional arguments passed to \code{ggplot2::ggsave}.
#' @return The paths of the written files, invisibly.
#' @export
SavePlotList <- function(plot_list, outdir = NULL, file = NULL, width = 6, height = 6, device = 'pdf', ...) {
  if (is.null(outdir) == is.null(file)) {
    stop('Specify exactly one of outdir or file.')
  }
  if (length(plot_list) == 0) {
    return(invisible(character(0)))
  }
  if (!is.null(file)) {
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    on.exit(grDevices::dev.off())
    for (p in plot_list) print(p)
    return(invisible(file))
  }
  if (is.null(names(plot_list)) || any(!nzchar(names(plot_list)))) {
    names(plot_list) <- paste0('plot_', seq_along(plot_list))
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  paths <- file.path(outdir, paste0(names(plot_list), '.', device))
  for (i in seq_along(plot_list)) {
    ggplot2::ggsave(paths[i], plot_list[[i]], width = width, height = height, device = device, ...)
  }
  invisible(paths)
}


# Resolve the mods argument of the network plots: 'all' means every module except
# grey, otherwise check that the requested modules exist.
select_modules <- function(modules, mods, wgcna_name) {
  if (all('all' %in% mods)) {
    return(setdiff(levels(modules$module), 'grey'))
  }
  missing_mods <- setdiff(mods, as.character(modules$module))
  if (length(missing_mods) > 0) {
    stop(paste0('Modules not found in ', wgcna_name, ': ', paste(missing_mods, collapse = ', ')))
  }
  mods
}

# Keep a proportion of the edges of each edge color, either at random or the
# strongest ones, and rescale edge weights to 0-1 within each color.
subsample_edges <- function(edge_df, edge_prop, sample_edges) {
  edge_df <- do.call(rbind, lapply(split(edge_df, edge_df$color), function(cur_df) {
    n_keep <- round(nrow(cur_df) * edge_prop)
    keep <- if (sample_edges) {
      sample.int(nrow(cur_df), n_keep)
    } else {
      order(cur_df$value, decreasing = TRUE)[seq_len(n_keep)]
    }
    cur_df <- cur_df[keep, , drop = FALSE]
    # a single edge (or all-equal weights) would give NaN from scale01()
    cur_df$value <- if (length(unique(cur_df$value)) > 1) scale01(cur_df$value) else rep(1, nrow(cur_df))
    cur_df
  }))
  rownames(edge_df) <- NULL
  edge_df
}

deprecated_outdir <- function(fn) {
  warning(
    'The outdir argument of ', fn, '() is deprecated: it now returns a list of plots. ',
    'Save them with SavePlotList(plots, outdir = ...).',
    call. = FALSE
  )
}


# Edge colours for network plots: the module colour when both genes of an edge are in the
# same module, otherwise 'grey90'. Vectorised; replaces a per-edge future_sapply() loop that
# exported the whole TOM and Seurat object to every worker.
edge_colors <- function(edge_df, modules) {
  col1 <- as.character(modules$color[match(as.character(edge_df$Var1), modules$gene_name)])
  col2 <- as.character(modules$color[match(as.character(edge_df$Var2), modules$gene_name)])
  ifelse(col1 == col2, col1, 'grey90')
}
