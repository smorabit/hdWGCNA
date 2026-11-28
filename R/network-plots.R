#' ModuleNetworkPlot
#'
#' Visualizes the top hub genes for selected modules as a circular network plot with one
#' inner circle and one outer circle.
#'
#' @param seurat_obj A Seurat object
#' @param n_inner Number of genes to plot on the inner circle
#' @param n_outer Number of genes to plot on the outer circle.
#' @param n_conns Number of gene-gene co-expression connections to plot, sorted by the strongest connections.
#' @param mods Names of the modules to plot. If mods = "all", all modules are plotted.
#' @param outdir The directory where the plots will be stored.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param plot_size A vector containing the width and height of the network plots. example: c(5,5)
#' @param edge.alpha value between 0 and 1 determining the alpha (transparency) scaling factor for the network edges
#' @param edge.width value determining the width of the network edges
#' @param vertex.label.cex vertex label font size
#' @param vertex.size vertex size
#' @param ... Currently unused.
#' @return Called for its side effect of writing one PDF per module to \code{outdir}; returns \code{NULL} invisibly.
#' @export
ModuleNetworkPlot <- function(
  seurat_obj,
  n_inner = 10,
  n_outer = 15,
  n_conns = 500,
  mods = 'all',
  outdir = 'ModuleNetworks',
  wgcna_name = NULL,
  plot_size = c(6, 6),
  edge.alpha = 0.25,
  edge.width = 1,
  vertex.label.cex = 1,
  vertex.size = 6, ...
) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get modules, MEs:
  modules <- GetModules(seurat_obj, wgcna_name)

  # using all modules?
  if (all('all' %in% mods)) {
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']
  } else {
    # check that the modules are present
    if (!all(mods %in% unique(as.character(modules$module)))) {
      stop(paste0('Some selected modules are not found in wgcna_name: ', wgcna_name))
    }
    modules <- modules %>% subset(module %in% mods)
  }

  # check if we have eigengene-based connectivities:
  if (!all(paste0('kME_', as.character(mods)) %in% colnames(modules))) {
    stop('Eigengene-based connectivity (kME) not found. Did you run ModuleEigengenes and ModuleConnectivity?')
  }

  # create output folder
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }

  # tell the user that the output is going to the output dir
  cat(paste0('Writing output files to ', outdir))

  # get TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # get hub genes:
  n_hubs <- n_inner + n_outer
  hub_df <- GetHubGenes(seurat_obj, n_hubs = n_hubs, wgcna_name = wgcna_name)

  # loop over modules
  for (cur_mod in mods) {
    cur_color <- modules %>%
      subset(module == cur_mod) %>%
      .$color %>%
      unique()
    cur_genes <- subset(hub_df, module == cur_mod) %>% .$gene_name
    n_genes <- length(cur_genes)

    # skip if there's too few genes:
    if (n_genes < (n_inner + 1)) {
      print(paste0('Skipping ', cur_mod, ', too few genes to plot.'))
      next
    }
    print(cur_mod)

    # name of column with current kME info
    cur_kME <- paste0('kME_', cur_mod)

    # Identify the columns in the TOM that correspond to these hub genes
    matchind <- match(cur_genes, colnames(TOM))
    reducedTOM <- TOM[matchind, matchind]
    orderind <- order(reducedTOM, decreasing = TRUE)

    if (n_conns > ncol(reducedTOM)**2) {
      cur_n_conns <- ncol(reducedTOM)**2
    } else {
      cur_n_conns <- n_conns
    }

    # only  keep top connections
    connections2keep <- orderind[1:cur_n_conns]
    connections2drop <- orderind[n_conns + seq_along(orderind)]
    reducedTOM[connections2drop] <- 0

    # scale between 0 and 1
    reducedTOM <- scale01(reducedTOM)

    # melt TOM into long format
    edge_df <- reducedTOM %>% reshape2::melt()
    edge_df$color_alpha <- alpha(cur_color, alpha = edge_df$value)

    # top 10 as center
    gA <- igraph::graph_from_adjacency_matrix(as.matrix(reducedTOM[1:n_inner, 1:n_inner]), mode = 'undirected', weighted = TRUE, diag = FALSE)
    gB <- igraph::graph_from_adjacency_matrix(as.matrix(reducedTOM[(n_inner + 1):n_genes, (n_inner + 1):n_genes]), mode = 'undirected', weighted = TRUE, diag = FALSE)
    layoutCircle <- rbind(igraph::layout_in_circle(gA) / 2, igraph::layout_in_circle(gB))

    g1 <- igraph::graph_from_data_frame(
      edge_df,
      directed = FALSE
    )

    pdf(paste0(outdir, '/', cur_mod, '.pdf'), width = plot_size[1], height = plot_size[2], useDingbats = FALSE)
    plot(g1,
      edge.color = adjustcolor(igraph::E(g1)$color_alpha, alpha.f = edge.alpha),
      edge.curved = 0,
      edge.width = edge.width,
      # edge.color=adjustcolor(cur_color, alpha.f=0.25),
      # edge.alpha=edge.alpha,
      vertex.color = cur_color,
      vertex.label = as.character(cur_genes),
      vertex.label.dist = 1.1,
      vertex.label.degree = -pi / 4,
      vertex.label.color = 'black',
      vertex.label.family = 'Helvetica',
      vertex.label.font = 3,
      vertex.label.cex = vertex.label.cex,
      vertex.frame.color = 'black',
      layout = jitter(layoutCircle),
      vertex.size = vertex.size,
      main = paste(cur_mod)
    )
    dev.off()
  }
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
#' @param return_graph logical determining whether we return the graph (TRUE) or plot the graph (FALSE)
#' @param edge.alpha Scaling factor for the edge opacity
#' @param vertex.label.cex The font size of the gene labels
#' @param hub.vertex.size The size of the hub gene nodes
#' @param other.vertex.size The size of the other gene nodes
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Additional arguments passed to both \code{igraph::layout_with_fr} and \code{plot.igraph}.
#' @return If \code{return_graph = TRUE}, an igraph object. Otherwise the network is drawn on the current graphics device and \code{NULL} is returned invisibly.
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
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get modules, MEs:
  MEs <- GetMEs(seurat_obj, wgcna_name = wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)

  # using all modules?
  if (all('all' %in% mods)) {
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']
  } else {
    # check that the modules are present
    if (!all(mods %in% unique(as.character(modules$module)))) {
      stop(paste0('Some selected modules are not found in wgcna_name: ', wgcna_name))
    }
    modules <- modules %>% subset(module %in% mods)
  }

  # get TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # get hub genes:
  hub_list <- lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur[, c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(n_hubs) %>%
      .$gene_name
  })
  names(hub_list) <- mods

  # sample the same number of genes in each module
  other_genes <- modules %>%
    subset(!(gene_name %in% unlist(hub_list))) %>%
    group_by(module) %>%
    sample_n(n_other, replace = TRUE) %>%
    .$gene_name %>%
    unique()

  # subset TOM by the selected genes:
  selected_genes <- c(unlist(hub_list), other_genes)
  selected_modules <- modules %>% subset(gene_name %in% selected_genes)
  subset_TOM <- TOM[selected_genes, selected_genes]

  # setup for network plot
  selected_modules$geneset <- ifelse(
    selected_modules$gene_name %in% other_genes, 'other', 'hub'
  )
  selected_modules$size <- ifelse(selected_modules$geneset == 'hub', hub.vertex.size, other.vertex.size)
  selected_modules$label <- ifelse(selected_modules$geneset == 'hub', as.character(selected_modules$gene_name), '')
  selected_modules$fontcolor <- ifelse(selected_modules$color == 'black', 'gray50', 'black')

  # make sure all nodes have at least one edge!!
  edge_cutoff <- min(sapply(seq_len(nrow(subset_TOM)), function(i) {
    max(subset_TOM[i, ])
  }))
  edge_df <- reshape2::melt(subset_TOM) %>% subset(value >= edge_cutoff)

  edge_df$color <- edge_colors(edge_df, modules)

  # subset edges:
  groups <- unique(edge_df$color)
  print(groups)
  if (sample_edges) {
    # randomly sample
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_sample <- sample(1:n_edges, round(n_edges * edge_prop))
      cur_df[cur_sample, ]
    }))
  } else {
    # get top strongest edges
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_df %>% dplyr::top_n(round(n_edges * edge_prop), wt = value)
    }))
  }

  edge_df <- temp

  # scale edge values between 0 and 1 for each module
  edge_df <- edge_df %>%
    group_by(color) %>%
    mutate(value = scale01(value))

  edge_df$color <- sapply(seq_len(nrow(edge_df)), function(i) {
    a <- edge_df$value[i]
    # if(edge_df$value[i] < 0.05){a=0.05}
    alpha(edge_df$color[i], alpha = a)
  })


  g <- igraph::graph_from_data_frame(
    edge_df,
    directed = FALSE,
    vertices = selected_modules
  )

  l <- igraph::layout_with_fr(g, ...)

  if (return_graph) {
    return(g)
  }

  plot(
    g,
    layout = l,
    edge.color = adjustcolor(igraph::E(g)$color, alpha.f = edge.alpha),
    vertex.size = igraph::V(g)$size,
    edge.curved = 0,
    edge.width = 0.5,
    vertex.color = igraph::V(g)$color,
    vertex.frame.color = igraph::V(g)$color,
    vertex.label = igraph::V(g)$label,
    vertex.label.family = 'Helvetica', # vertex.label.font=vertex_df$font,
    vertex.label.font = 3,
    vertex.label.color = igraph::V(g)$fontcolor,
    vertex.label.cex = vertex.label.cex,
    ...
  )
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
#' Makes a igraph network plot using the module UMAP
#'
#' @param seurat_obj A Seurat object
#' @param sample_edges logical determining whether we downsample edges for plotting (TRUE), or take the strongst edges.
#' @param edge_prop proportion of edges to plot. If sample_edges=FALSE, the strongest edges are selected.
#' @param label_hubs the number of hub genes to label in each module
#' @param edge.alpha scaling factor for edge opacity
#' @param vertex.label.cex font size for labeled genes
#' @param label_genes Optional character vector of additional genes to label, alongside the top hub genes.
#' @param return_graph logical determining whether to plot thr graph (FALSE) or return the igraph object (TRUE)
#' @param keep_grey_edges logical determining whether to show edges between genes in different modules (grey edges)
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Currently unused.
#' @return If \code{return_graph = TRUE}, an igraph object. Otherwise the network is drawn on the current graphics device and \code{NULL} is returned invisibly.
#' @keywords scRNA-seq
#' @export
ModuleUMAPPlot <- function(
  seurat_obj,
  sample_edges = TRUE, # TRUE if we sample edges randomly, FALSE if we take the top edges
  edge_prop = 0.2,
  label_hubs = 5, # how many hub genes to label?
  edge.alpha = 0.25,
  vertex.label.cex = 0.5,
  label_genes = NULL,
  return_graph = FALSE, # this returns the igraph object instead of plotting
  keep_grey_edges = TRUE,
  wgcna_name = NULL,
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get the TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # get modules,
  modules <- GetModules(seurat_obj, wgcna_name)

  # get the UMAP df:
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  mods <- levels(umap_df$module)
  mods <- mods[mods != 'grey']

  # subset the TOM:
  subset_TOM <- TOM[umap_df$gene, umap_df$gene[umap_df$hub == 'hub']]

  # genes to label:
  # hub_labels <- selected_modules %>% group_by(module) %>% top_n(label_hubs, wt=kME) %>% .$gene_name
  hub_list <- lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur[, c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(label_hubs) %>%
      .$gene_name
  })
  names(hub_list) <- mods
  hub_labels <- as.character(unlist(hub_list))
  print('hub labels')
  print(hub_labels)
  print(label_genes)
  if (is.null(label_genes)) {
    label_genes <- hub_labels
  } else {
    if (!any(label_genes %in% umap_df$gene)) {
      stop('Some genes in label_genes not found in the UMAP.')
    }
    label_genes <- unique(c(label_genes, hub_labels))
  }
  print(label_genes)

  # subset module df by genes in the UMAP df:
  selected_modules <- modules[umap_df$gene, ]
  selected_modules <- cbind(selected_modules, umap_df[, c('UMAP1', 'UMAP2', 'hub', 'kME')])

  selected_modules$label <- ifelse(selected_modules$gene_name %in% label_genes, selected_modules$gene_name, '')
  selected_modules$fontcolor <- ifelse(selected_modules$color == 'black', 'gray50', 'black')

  # set frome color
  # same color as module for all genes, black outline for the selected hub genes
  selected_modules$framecolor <- ifelse(selected_modules$gene_name %in% label_genes, 'black', selected_modules$color)

  # melt TOM into long format
  edge_df <- subset_TOM %>% reshape2::melt()
  print(dim(edge_df))

  # set color of each edge based on value:
  edge_df$color <- edge_colors(edge_df, selected_modules)

  # keep grey edges?
  if (!keep_grey_edges) {
    edge_df <- edge_df %>% subset(color != 'grey90')
  }

  # subset edges:
  groups <- unique(edge_df$color)
  if (sample_edges) {
    # randomly sample
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_sample <- sample(1:n_edges, round(n_edges * edge_prop))
      cur_df[cur_sample, ]
    }))
  } else {
    # get top strongest edges
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_df %>% dplyr::top_n(round(n_edges * edge_prop), wt = value)
    }))
  }

  edge_df <- temp
  print(dim(edge_df))

  # scale edge values between 0 and 1 for each module
  edge_df <- edge_df %>%
    group_by(color) %>%
    mutate(value = scale01(value))

  # edges & vertices are plotted in igraph starting with the first row, so re-order s.t. strong edges are on bottom, all gray on the top of the table:
  edge_df <- edge_df %>% arrange(value)
  edge_df <- rbind(
    subset(edge_df, color == 'grey90'),
    subset(edge_df, color != 'grey90')
  )

  # set alpha of edges based on kME
  edge_df$color_alpha <- ifelse(
    edge_df$color == 'grey90',
    alpha(edge_df$color, alpha = edge_df$value / 2),
    alpha(edge_df$color, alpha = edge_df$value)
  )

  # re-order vertices so hubs are plotted on top
  selected_modules <- rbind(
    subset(selected_modules, hub == 'other'),
    subset(selected_modules, hub != 'other')
  )

  # re-order vertices so labeled genes are on top
  selected_modules <- rbind(
    subset(selected_modules, label == ''),
    subset(selected_modules, label != '')
  )

  # setup igraph:
  g <- igraph::graph_from_data_frame(
    edge_df,
    directed = FALSE,
    vertices = selected_modules
  )

  # print('making net')
  # print(head(edge_df))
  # print(head(selected_modules))

  if (return_graph) {
    return(g)
  }

  plot(
    g,
    layout = as.matrix(selected_modules[, c('UMAP1', 'UMAP2')]),
    # edge.color=adjustcolor(igraph::E(g)$color, alpha.f=edge.alpha),
    edge.color = adjustcolor(igraph::E(g)$color_alpha, alpha.f = edge.alpha),
    vertex.size = igraph::V(g)$kME * 3,
    edge.curved = 0,
    edge.width = 0.5,
    vertex.color = igraph::V(g)$color,
    vertex.label = igraph::V(g)$label,
    vertex.label.dist = 1.1,
    vertex.label.degree = -pi / 4,
    vertex.label.family = 'Helvetica', # vertex.label.font=vertex_df$font,
    vertex.label.font = 3,
    vertex.label.color = igraph::V(g)$fontcolor,
    vertex.label.cex = vertex.label.cex,
    vertex.frame.color = igraph::V(g)$framecolor,
    margin = 0
  )
}


#' ModuleCorrNetworks
#'
#' Plot Module Eigengene correlogram
#'
#' @param seurat_obj A Seurat object
#' @param cluster_col Column in \code{seurat_obj@meta.data} with cluster labels used to position modules. Defaults to \code{Idents(seurat_obj)}.
#' @param exclude_grey logical determining whether to exclude the grey module.
#' @param features Which module features to correlate: 'hMEs', 'MEs', 'scores', or 'average'.
#' @param reduction Name of the dimensionality reduction used to compute cluster centroids. Default = 'umap'.
#' @param cor_cutoff Minimum absolute correlation for an edge to be drawn. Default = 0.2.
#' @param label_vertices logical determining whether to label the module nodes.
#' @param edge_scale Scaling factor applied to the absolute correlation to set edge width.
#' @param vertex_size Size of the module nodes.
#' @param niter Number of iterations for the Fruchterman-Reingold layout.
#' @param vertex_frame logical determining whether to draw a black outline around the nodes (TRUE) or use the module color (FALSE).
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return Draws the module correlation network on the current graphics device and returns \code{NULL} invisibly.
#' @keywords scRNA-seq
#' @export
ModuleCorrNetwork <- function(
  seurat_obj, cluster_col = NULL, exclude_grey = TRUE,
  features = 'hMEs',
  reduction = 'umap', cor_cutoff = 0.2, label_vertices = FALSE, edge_scale = 5,
  vertex_size = 15, niter = 100, vertex_frame = FALSE, wgcna_name = NULL
) {
  check_suggests('qgraph', feature = 'ModuleCorrNetwork')

  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # Get module eigengenes
  if (features == 'hMEs') {
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  } else if (features == 'MEs') {
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  } else if (features == 'scores') {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else if (features == 'average') {
    MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
    restrict_range <- FALSE
  } else {
    (
      stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
    )
  }

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)

  # exclude grey
  if (exclude_grey) {
    MEs <- MEs[, colnames(MEs) != 'grey']
    modules <- modules %>% subset(color != 'grey')
  }

  # get list of modules
  mods <- colnames(MEs)

  # what clusters are we using?
  if (is.null(cluster_col)) {
    clusters <- Idents(seurat_obj)
  } else {
    clusters <- droplevels(seurat_obj@meta.data[[cluster_col]])
  }

  MEs$cluster <- clusters

  # compute average MEs for each cluster
  cluster_ME_av <- do.call(
    rbind, lapply(
      split(MEs, MEs$cluster),
      function(x) {
        colMeans(x[, mods])
      }
    )
  ) %>% as.data.frame()

  # which cluster mas highest expression of each module?
  top_clusters <- sapply(mods, function(x) {
    rownames(cluster_ME_av[cluster_ME_av[, x] == max(cluster_ME_av[, x]), ])
  })

  # get UMAP / tSNE centroids for these clusters to use as starting coordinates
  check_reduction(seurat_obj, reduction)
  red_df <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
  red_df$cluster <- clusters

  # compute average coords for each cluster
  red_av <- do.call(
    rbind, lapply(
      split(red_df, red_df$cluster),
      function(x) {
        colMeans(x[, 1:2])
      }
    )
  ) %>% as.data.frame()


  # get correlation matrix
  cor_mat <- Hmisc::rcorr(as.matrix(MEs[, seq_len(ncol(MEs)) - 1]))$r
  cor_mat[lower.tri(cor_mat)] <- NA

  # melt matrix and remove NAs
  cor_df <- reshape2::melt(cor_mat) %>% na.omit()

  # remove self-edges
  cor_df <- cor_df %>% subset(!(Var1 == Var2))

  # remove weak edges:
  cor_df <- cor_df %>% subset(abs(value) >= cor_cutoff)

  # vertex df:
  v_df <- data.frame(
    name = mods,
    cluster = as.character(top_clusters)
  )

  # add module colors:
  unique_mods <- dplyr::distinct(modules[, c('module', 'color')])
  rownames(unique_mods) <- unique_mods$module
  v_df$color <- unique_mods[v_df$name, 'color']

  # add reduction coords:
  v_df$x <- red_av[v_df$cluster, 1]
  v_df$y <- red_av[v_df$cluster, 2]

  # make graph:
  g <- igraph::graph_from_data_frame(
    cor_df,
    directed = FALSE,
    vertices = v_df
  )

  # igraph layout
  e <- igraph::as_edgelist(g, names = FALSE)
  # l <- igraph::layout_with_fr(
  #   g,
  #   weights=igraph::E(g)$value,
  #   coords = as.matrix(data.frame(x=igraph::V(g)$x, y=igraph::V(g)$y)),
  #   niter=niter
  # )
  #
  l <- qgraph::qgraph.layout.fruchtermanreingold(
    e,
    vcount = vcount(g),
    weights = igraph::E(g)$value,
    repulse.rad = (vcount(g)),
    # cool.exp = 0.5,
    # init = as.matrix(data.frame(x=igraph::V(g)$x, y=igraph::V(g)$y)),
    niter = niter,
    # max.delta = vcount(g)/2
  )


  # ggplot just to get the colors
  plot_df <- rbind(cor_df, data.frame(Var1 = c('x', 'y'), Var2 = c('y', 'x'), value = c(-1, 1)))
  temp <- ggplot(plot_df, aes(x = value, y = value, color = value)) +
    geom_point() +
    scale_color_gradient2(high = 'darkorchid1', mid = 'white', low = 'seagreen', midpoint = 0)
  temp <- ggplot_build(temp)
  igraph::E(g)$color <- temp$data[[1]]$colour[seq_len(nrow(cor_df))]

  # label the vertices?
  if (label_vertices) {
    labels <- igraph::V(g)$name
  } else {
    labels <- NA
  }

  # vertex_frame
  if (vertex_frame) {
    frame_color <- 'black'
  } else {
    frame_color <- igraph::V(g)$color
  }

  # plot the graph
  plot(
    g,
    layout = l,
    edge.color = igraph::E(g)$color,
    edge.curved = 0,
    edge.width = abs(igraph::E(g)$value) * edge_scale,
    vertex.color = igraph::V(g)$color,
    vertex.frame.color = frame_color,
    vertex.label = labels,
    vertex.label.family = 'Helvetica',
    vertex.label.color = 'black',
    vertex.label.cex = 0.5,
    vertex.size = vertex_size
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
