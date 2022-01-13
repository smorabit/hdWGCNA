
#' PlotDendrogram
#'
#' Plot WGCNA dendrogram
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' PlotDendrogram
PlotDendrogram <- function(
  seurat_obj, groupLabels="Module colors", wgcna_name=NULL,
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = "", ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get WGCNA network and module data
  net <- GetNetworkData(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)

  # plot dendrogram
  WGCNA::plotDendroAndColors(
    net$dendrograms[[1]],
    as.character(modules$color),
    groupLabels=groupLabels,
    dendroLabels = dendroLabels,
    hang = hang,
    addGuide = addGuide,
    guideHang = guideHang,
    main = main,
    ...
  )
}

#' ModuleCorrelogram
#'
#' Plot Module Eigengene correlogram
#'
#' @param seurat_obj A Seurat object
#' @param features What to plot? Can select hMEs, MEs, scores, or average
#' @keywords scRNA-seq
#' @export
#' @examples
#' MECorrelogram
ModuleCorrelogram <- function(
  seurat_obj, wgcna_name = NULL, exclude_grey=TRUE,
  features = 'hMEs',
  type = 'upper', order='original', method='ellipse',
  tl.col = 'black', tl.srt=45,
  sig.level = c(0.0001, 0.001, 0.01, 0.05),
#  insig='label_sig',
  pch.cex=0.7, col=NULL, ncolors=200, ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get MEs, module data from seurat object
  if(features == 'hMEs'){
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  } else if(features == 'MEs'){
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  } else if(features == 'scores'){
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else if(features == 'average'){
    MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
    restrict_range <- FALSE
  } else(
    stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
  )

  # convert to matrix
  MEs <- as.matrix(MEs)

  # exclude grey
  if(exclude_grey){
    MEs <- MEs[,colnames(MEs) != 'grey']
  }

  # setup color scheme:
  if(is.null(col)){
    colfunc <- grDevices::colorRampPalette(c('seagreen',  'white', 'darkorchid1'))
    col = colfunc(ncolors)
  }

  # perform correlation
  res <- Hmisc::rcorr(MEs)
  resP <- corrplot::cor.mtest(MEs, conf.level=0.95)$p

  # plot correlogram
  corrplot::corrplot(
    res$r, p.mat = resP,
    type=type, order=order,
    method=method, tl.col=tl.col,
    tl.srt=tl.srt, sig.level=sig.level,
    # insig=insig,
    pch.cex=pch.cex,
    col = col, ...
  )

}


#' ModuleCorrNetworks
#'
#' Plot Module Eigengene correlogram
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleCorrNetwork
ModuleCorrNetwork <- function(
  seurat_obj, wgcna_name=NULL, cluster_col=NULL, exclude_grey=TRUE,
  features = 'hMEs',
  reduction='umap', cor_cutoff=0.2, label_vertices=FALSE, edge_scale=5,
  vertex_size=15, niter=100, vertex_frame=FALSE
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # Get module eigengenes
  if(features == 'hMEs'){
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  } else if(features == 'MEs'){
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  } else if(features == 'scores'){
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else if(features == 'average'){
    MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
    restrict_range <- FALSE
  } else(
    stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
  )

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)

  # exclude grey
  if(exclude_grey){
    MEs <- MEs[,colnames(MEs) != 'grey']
    modules <- modules %>% subset(color != 'grey')
  }

  # get list of modules
  mods  <-  colnames(MEs)

  # what clusters are we using?
  if(is.null(cluster_col)){
    clusters <- Idents(seurat_obj)
  } else{
    clusters <- droplevels(seurat_obj@meta.data[[cluster_col]])
  }

  MEs$cluster <- clusters

  # compute average MEs for each cluster
  cluster_ME_av <- do.call(
    rbind, lapply(
      split(MEs, MEs$cluster),
      function(x){colMeans(x[,mods])
  })) %>% as.data.frame

  # which cluster mas highest expression of each module?
  top_clusters <- sapply(mods, function(x){
    rownames(cluster_ME_av[cluster_ME_av[,x] == max(cluster_ME_av[,x]),])
  })

  # get UMAP / tSNE centroids for these clusters to use as starting coordinates
  red_df <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
  red_df$cluster <- clusters

  # compute average coords for each cluster
  red_av <- do.call(
    rbind, lapply(
      split(red_df, red_df$cluster),
      function(x){colMeans(x[,1:2])
  })) %>% as.data.frame



  # get correlation matrix
  cor_mat <- Hmisc::rcorr(as.matrix(MEs[,1:ncol(MEs)-1]))$r
  cor_mat[lower.tri(cor_mat)] <- NA

  # melt matrix and remove NAs
  cor_df <- reshape2::melt(cor_mat) %>% na.omit

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
  unique_mods <- distinct(modules[,c('module', 'color')])
  rownames(unique_mods) <- unique_mods$module
  v_df$color <- unique_mods[v_df$name,'color']

  # add reduction coords:
  v_df$x <- red_av[v_df$cluster, 1]
  v_df$y <- red_av[v_df$cluster, 2]

  # make graph:
  g <- igraph::graph_from_data_frame(
    cor_df,
    directed=FALSE,
    vertices=v_df
  )

  # igraph layout
  e <-  get.edgelist(g, name=FALSE)
  # l <- igraph::layout_with_fr(
  #   g,
  #   weights=E(g)$value,
  #   coords = as.matrix(data.frame(x=V(g)$x, y=V(g)$y)),
  #   niter=niter
  # )
  #
  l <- qgraph::qgraph.layout.fruchtermanreingold(
    e, vcount = vcount(g),
    weights=E(g)$value,
    repulse.rad=(vcount(g)),
    #cool.exp = 0.5,
    # init = as.matrix(data.frame(x=V(g)$x, y=V(g)$y)),
    niter=niter,
    #max.delta = vcount(g)/2
  )


  # ggplot just to get the colors
  plot_df <- rbind(cor_df, data.frame(Var1=c('x', 'y'), Var2=c('y', 'x'), value=c(-1,1)))
  temp <- ggplot(plot_df, aes(x=value, y=value, color=value)) +
    geom_point() + scale_color_gradient2(high='darkorchid1', mid='white', low='seagreen', midpoint=0)
  temp <- ggplot_build(temp)
  E(g)$color <- temp$data[[1]]$colour[1:nrow(cor_df)]

  # label the vertices?
  if(label_vertices){labels <- V(g)$name} else{labels <- NA}

  # vertex_frame
  if(vertex_frame){frame_color <- 'black'} else{frame_color <- V(g)$color}

  # plot the graph
  plot(
    g, layout=l,
    edge.color=E(g)$color,
    edge.curved=0,
    edge.width=abs(E(g)$value) * edge_scale,
    vertex.color=V(g)$color,
    vertex.frame.color=frame_color,
    vertex.label=labels,
    vertex.label.family='Helvetica',
    vertex.label.color = 'black',
    vertex.label.cex=0.5,
    vertex.size=vertex_size
  )
}


#' ModuleFeaturePlot
#'
#' Plot module eigengenes as a FeaturePlot
#'
#' @param seurat_obj A Seurat object
#' @param features What to plot? Can select hMEs, MEs, scores, or average
#' @param order TRUE, FALSE, or "shuffle" are valid options
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleFeaturePlot
ModuleFeaturePlot<- function(
  seurat_obj, module_names=NULL, wgcna_name = NULL,
  reduction='umap', features = 'hMEs',
  order=TRUE, restrict_range=TRUE, point_size = 0.5, alpha=1,
  label_legend = FALSE, ucell = FALSE
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get MEs, module data from seurat object
  if(features == 'hMEs'){
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  } else if(features == 'MEs'){
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  } else if(features == 'scores'){
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else if(features == 'average'){
    MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
    restrict_range <- FALSE
  } else(
    stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
  )

  # override restrict_range if ucell
  if(ucell){restrict_range <- FALSE}

  # use all modules except gray if not specified by the user
  modules <- GetModules(seurat_obj, wgcna_name)
  if(is.null(module_names)){
    module_names <- colnames(MEs)
    module_names <- module_names[module_names != 'grey']
  }

  # get  reduction from seurat obj
  umap <- seurat_obj@reductions[[reduction]]@cell.embeddings
  x_name <- colnames(umap)[1]
  y_name <- colnames(umap)[2]

  # merge into one df for plotting
  plot_df <- cbind(umap, MEs) %>% as.data.frame()

  plot_list <- list()
  for(cur_mod in module_names){

    # get the color for this module:
    cur_color <- modules %>% subset(module == cur_mod) %>% .$color %>% unique

    # reset the range of the plot:
    plot_range <- plot_df[,cur_mod] %>% range
    if(restrict_range){
      if(abs(plot_range[1]) > abs(plot_range[2])){
        plot_range[1] <- -1*plot_range[2]
      } else{
        plot_range[2] <- -1*plot_range[1]
      }
      plot_df[,cur_mod] <- ifelse(plot_df[,cur_mod] > plot_range[2], plot_range[2], plot_df[,cur_mod])
      plot_df[,cur_mod] <- ifelse(plot_df[,cur_mod] < plot_range[1], plot_range[1], plot_df[,cur_mod])
    }

    # order points:
    if(order == TRUE){
      plot_df <- plot_df %>% dplyr::arrange_(cur_mod)
    } else if(order == "shuffle"){
      plot_df <- plot_df[sample(nrow(plot_df)),]
    }

    # label for plot:
    label <- cur_mod

    # plot with ggplot
     p <- plot_df %>%
      ggplot(aes_string(x=x_name, y=y_name, color=cur_mod)) +
      geom_point(size=point_size, alpha=alpha) +
      ggtitle(label) + umap_theme +
      labs(color="")

    # UCell?
    if(!ucell){
      p <- p + scale_color_gradient2(
        low='grey75', mid='grey95', high=cur_color,
        breaks = c(plot_range[1], plot_range[2]),
        labels = c('-', '+'),
        guide = guide_colorbar(ticks=FALSE, barwidth=0.5, barheight=4)
      )
    } else{
      p <- p + scale_color_gradient(
        low='grey95', high=cur_color,
        breaks = c(plot_range[1], plot_range[2]),
        labels = c('0', '+'),
        guide = guide_colorbar(ticks=FALSE, barwidth=0.5, barheight=4)
      )
    }

    plot_list[[cur_mod]] <- p

  }

  # return plot
  if(length(plot_list) == 1){
    p <- plot_list[[1]]
  } else{
    p <- plot_list
  }

  p

}





#' EnrichrBarPlot
#'
#' Makes barplots from Enrichr data
#'
#' @param seurat_obj A Seurat object
#' @param dbs List of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' EnrichrBarPlot
EnrichrBarPlot <- function(
  seurat_obj, outdir = "enrichr_plots",
  n_terms = 25, plot_size = c(6,15),
  wgcna_name=NULL, logscale=FALSE, ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get Enrichr table
  enrichr_df <- GetEnrichrTable(seurat_obj, wgcna_name)

  # helper function to wrap text
  wrapText <- function(x, len) {
      sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
  }

  # make output dir if it doesn't exist:
  if(!dir.exists(outdir)){dir.create(outdir)}

  # loop through modules:
  for(i in 1:length(mods)){

    cur_mod <- mods[i]
    cur_terms <- subset(enrichr_df, module == cur_mod)
    print(cur_mod)

    # get color for this module:
    cur_color <- modules %>% subset(module == cur_mod) %>% .$color %>% unique %>% as.character

    # skip if there are not any terms for this module:
    if(nrow(cur_terms) == 0){next}
    cur_terms$wrap <- wrapText(cur_terms$Term, 45)

    # plot top n_terms as barplot
    plot_list <- list()
    for(cur_db in dbs){

      plot_df <- subset(cur_terms, db==cur_db) %>% top_n(n_terms, wt=Combined.Score)

      # text color:
      if(cur_mod == 'black'){
        text_color = 'grey'
      } else {
        text_color = 'black'
      }

      # logscale?
      if(logscale){
        plot_df$Combined.Score <- log(plot_df$Combined.Score)
        lab <- 'Enrichment log(combined score)'
        x <- 0.2
      } else{lab <- 'Enrichment (combined score)'; x <- 5}

      # make bar plot:
      plot_list[[cur_db]] <- ggplot(plot_df, aes(x=Combined.Score, y=reorder(wrap, Combined.Score)))+
        geom_bar(stat='identity', position='identity', color='white', fill=cur_color) +
        geom_text(aes(label=wrap), x=x, color=text_color, size=3.5, hjust='left') +
        ylab('Term') + xlab(lab) + ggtitle(cur_db) +
        theme(
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.title = element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          plot.title = element_text(hjust = 0.5)
        )
    }

    # make pdfs in output dir
    pdf(paste0(outdir, '/', cur_mod, '.pdf'), width=plot_size[1], height=plot_size[2])
    for(plot in plot_list){
      print(plot)
    }
    dev.off()
  }
}


#' EnrichrDotPlot
#'
#' Makes barplots from Enrichr data
#'
#' @param seurat_obj A Seurat object
#' @param dbs List of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' EnrichrBarPlot
EnrichrDotPlot <- function(
  seurat_obj, database, mods="all", outdir = "enrichr_plots",
  n_terms = 3, break_ties=TRUE,
  wgcna_name=NULL, logscale=TRUE, ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)

  # using all modules?
  if(mods == 'all'){
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']
  }

  # get Enrichr table
  enrichr_df <- GetEnrichrTable(seurat_obj, wgcna_name)

  # helper function to wrap text
  wrapText <- function(x, len) {
      sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
  }

  # get data to plot
  plot_df <- enrichr_df %>%
    subset(db == database & module %in% mods) %>%
    group_by(module) %>%
    top_n(n_terms, wt=Combined.Score)

  print(table(plot_df$db))

  # sometimes top_n returns more than the desired number if there are ties. so here
  # we just randomly sample to break ties:
  if(break_ties){
    plot_df <- do.call(rbind, lapply(plot_df %>% group_by(module) %>% group_split, function(x){x[sample(n_terms),]}))
  }

  plot_df$Term <- wrapText(plot_df$Term, 45)

  # set modules factor and re-order:
  plot_df$module <- factor(
    as.character(plot_df$module),
    levels=levels(modules$module)
  )
  plot_df <- arrange(plot_df, module)

  # set Terms factor
  plot_df$Term <- factor(
    as.character(plot_df$Term),
    levels = unique(as.character(plot_df$Term))
  )

  # logscale?
  if(logscale){
    plot_df$Combined.Score <- log(plot_df$Combined.Score)
    lab <- 'Enrichment\nlog(combined score)'
    x <- 0.2
  } else{lab <- 'Enrichment\n(combined score)'; x <- 5}


  p <- plot_df  %>%
    ggplot(aes(x=module, y=rev(Term))) +
    geom_point(aes(size=Combined.Score), color=plot_df$module) +
    RotatedAxis() +
    ylab('') + xlab('') + labs(size=lab) +
    ggtitle(database) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )

  p
}


#' ModuleNetworkPlot
#'
#' Makes barplots from Enrichr data
#'
#' @param seurat_obj A Seurat object
#' @param dbs List of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleNetworkPlot
ModuleNetworkPlot <- function(
  seurat_obj, mods="all", outdir="ModuleNetworks",
  plot_size = c(6,6), wgcna_name=NULL,
  edge.alpha=0.25, vertex.label.cex=1, vertex.size=6, ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules, MEs:
  MEs <- GetMEs(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)

  # using all modules?
  if(mods == 'all'){
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']
  }

  # create output folder
  if(!dir.exists(outdir)){dir.create(outdir)}

  # get TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # loop over modules
  for(cur_mod in mods){

    cur_color <- modules %>% subset(module == cur_mod) %>% .$color %>% unique

    # number of genes, connections
    # might make a setting to change this later but I will also have to change
    # how the graph layout works
    n_genes = 25;
    n_conns = 500;

    # name of column with current kME info
    cur_kME <- paste0('kME_', cur_mod)

    # get genes in this module:
    cur <- subset(modules, module == cur_mod)
    rowind <- cur[,c('gene_name', cur_kME)] %>%
      top_n(n_genes) %>% .$gene_name
    cur <- cur[rowind,]

    # arrange by cur_kME
    cur <- cur %>% arrange(get(cur_kME), descending=TRUE)

    if (nrow(cur) < n_genes) {
      n_genes <-  nrow(cur);
      n_conns <- n_genes * (n_genes - 1);
    }

    # Identify the columns in the TOM that correspond to these hub genes
    matchind <- match(cur$gene_name, colnames(TOM))
    reducedTOM = TOM[matchind,matchind]
    orderind <- order(reducedTOM,decreasing=TRUE)

    # only  keep top connections
    connections2keep <- orderind[1:n_conns];
    reducedTOM <- matrix(0,nrow(reducedTOM),ncol(reducedTOM));
    reducedTOM[connections2keep] <- 1;

    # top 10 as center
    gA <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
    gB <- graph.adjacency(as.matrix(reducedTOM[11:n_genes,11:n_genes]),mode="undirected",weighted=TRUE,diag=FALSE)
    layoutCircle <- rbind(layout.circle(gA)/2,layout.circle(gB))

    g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)

    pdf(paste0(outdir, '/', cur_mod,'.pdf'), width=plot_size[1], height=plot_size[2], useDingbats=FALSE);
    plot(g1,
      edge.color=adjustcolor(cur_color, alpha.f=0.25),
      edge.alpha=edge.alpha,
      vertex.color=cur_color,
      vertex.label=as.character(cur$gene_name),
      vertex.label.dist=1.1,
      vertex.label.degree=-pi/4,
      vertex.label.color="black",
      vertex.label.family='Helvetica',
      vertex.label.font = 3,
      vertex.label.cex=vertex.label.cex,
      vertex.frame.color='black',
      layout= jitter(layoutCircle),
      vertex.size=vertex.size,
      main=paste(cur_mod)
    )
    dev.off();

  }

}



#' HubGeneNetworkPlot
#'
#' Makes barplots from Enrichr data
#'
#' @param seurat_obj A Seurat object
#' @param dbs List of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' HubGeneNetworkPlot
HubGeneNetworkPlot <- function(
  seurat_obj, mods="all", n_hubs=6, n_other=3,
  plot_size = c(6,6), wgcna_name=NULL,
  edge.alpha=0.25, vertex.label.cex=0.5, hub.vertex.size=4,
  other.vertex.size=1, repulse.exp=3,  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules, MEs:
  MEs <- GetMEs(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)

  # using all modules?
  if(mods == 'all'){
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']
  } else{
    modules <- modules %>% subset(module %in% mods)
  }

  # get TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # get hub genes:
  hub_list <- lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(n_hubs) %>% .$gene_name
  })
  names(hub_list) <- mods

  # sample the same number of genes in each module
  other_genes <- modules %>%
    subset(!(gene_list %in% unlist(hub_list))) %>%
    group_by(module) %>%
    sample_n(n_other, replace=TRUE) %>%
    .$gene_name %>% unique

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
  edge_cutoff <- min(sapply(1:nrow(subset_TOM), function(i){max(subset_TOM[i,])}))
  edge_df <- subset_TOM %>% melt %>% subset(value >= edge_cutoff)

  # remove nodes with fewer than n edges:
  n_fewer = 5
  remove_nodes <- table(edge_df$Var1)[table(edge_df$Var1) < n_fewer] %>% names
  edge_df <- subset(edge_df, !(Var1 %in% remove_nodes) & !(Var2 %in% remove_nodes))
  selected_modules <- subset(selected_modules, !(gene_name %in% remove_nodes))

  # scale edge values between 0 and 1
  edge_df$value <- (edge_df$value - min(edge_df$value)) / (max(edge_df$value) - min(edge_df$value))

  # set color of each edge based on value:
  edge_df$color <- sapply(1:nrow(edge_df), function(i){
    gene1 = as.character(edge_df[i,'Var1'])
    gene2 = as.character(edge_df[i,'Var2'])

    col1 <- modules %>% subset(gene_name == gene1) %>% .$color
    col2 <- modules %>% subset(gene_name == gene2) %>% .$color

    if(col1 == col2){
      col = col1
    } else{
      col = 'gray90'
    }
    col
  })

  edge_df$color <- sapply(1:nrow(edge_df), function(i){
    a = edge_df$value[i]
    #if(edge_df$value[i] < 0.05){a=0.05}
    alpha(edge_df$color[i], alpha=a)
  })

  g <- igraph::graph_from_data_frame(
    edge_df,
    directed=FALSE,
    vertices=selected_modules
  )

  # qgraph layout
  e <-  get.edgelist(g, name=FALSE)
  l <- qgraph.layout.fruchtermanreingold(
    e, vcount = vcount(g),
    weights=edge_df$value,
    repulse.rad=(vcount(g)^repulse.exp),
    #cool.exp = 0.5,
    niter=2500,
    #max.delta = vcount(g)/2
  )

  # make a communities? (nope, this looks bad)
  # comm <- igraph::make_clusters(
  #   g,
  #   membership=as.numeric(as.factor(V(g)$module)),
  #   algorithm="scWGCNA"
  # )

  plot(
    g, layout=l,
    edge.color=adjustcolor(E(g)$color, alpha.f=edge.alpha),
    vertex.size=V(g)$size,
    edge.curved=0,
    edge.width=0.5,
    vertex.color=V(g)$color,
    vertex.frame.color=V(g)$color,
    vertex.label=V(g)$label,
    vertex.label.family='Helvetica', #vertex.label.font=vertex_df$font,
    vertex.label.font = 3,
    vertex.label.color = V(g)$fontcolor,
    vertex.label.cex=vertex.label.cex,
    ...
  )

}


#' OverlapDotPlot
#'
#' Makes barplots from Enrichr data
#'
#' @param seurat_obj A Seurat object
#' @param dbs List of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' OverlapDotPlot
OverlapDotPlot <- function(
  overlap_df, plot_var = 'odds_ratio',
  logscale=TRUE, neglog=FALSE, plot_significance=TRUE, ...
){

  label <- plot_var
  if(logscale){
    overlap_df[[plot_var]] <- log(overlap_df[[plot_var]])
    label <- paste0('log(', plot_var, ')')
  }
  if(neglog){
    overlap_df[[plot_var]] <- -1 * overlap_df[[plot_var]]
    label <- paste0('-', label)
  }

  p <- overlap_df %>% ggplot(aes(x=module, y=group)) +
    geom_point(aes(
      size=get(plot_var),
      alpha=get(plot_var)),
      color=overlap_df$color
    ) +
    RotatedAxis() +
    ylab('') + xlab('') + labs(size=label, alpha=label) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )

  # plot significance level?
  if(plot_significance){
    p <- p + geom_text(aes(label=Significance))
  }

  p
}

#' OverlapBarPlot
#'
#' Makes barplots from Enrichr data
#'
#' @param seurat_obj A Seurat object
#' @param dbs List of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' OverlapBarPlot
OverlapBarPlot <- function(
  overlap_df, plot_var = 'odds_ratio',
  logscale=TRUE, neglog=FALSE, ...
){

  label <- plot_var

  if(plot_var == 'odds_ratio'){
    yint <- 1
  } else if(plot_var == 'fdr'){
    yint <- 0.05
  }

  if(logscale){
    overlap_df[[plot_var]] <- log(overlap_df[[plot_var]])
    label <- paste0('log(', plot_var, ')')
    yint = log(yint)
  }
  if(neglog){
    overlap_df[[plot_var]] <- -1 * overlap_df[[plot_var]]
    label <- paste0('-', label)
    yint = -1 * yint
  }

  groups <- overlap_df$group %>% as.character %>% unique

  plot_list <- list()
  for(cur_group in groups){

    cur_df <- overlap_df %>%
      subset(group == cur_group)

    p <- cur_df %>%
      ggplot(aes(x=reorder(module, get(plot_var)), y=get(plot_var))) +
      geom_bar(stat='identity', fill=cur_df$color) +
      coord_flip() +
      xlab('') + ylab(label) +
      ggtitle(cur_group) +
      theme(
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )

      if(plot_var == 'fdr' | plot_var == 'odds_ratio'){
        p <- p + geom_hline(yintercept=yint, linetype='dashed', color='gray')
      }

      plot_list[[cur_group]] <- p
  }

  plot_list

}
