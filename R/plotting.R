
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
  seurat_obj, MEs2=NULL,
  features = 'hMEs',
  order='original', method='ellipse',
  exclude_grey=TRUE, type='upper',
  tl.col = 'black', tl.srt=45,
  sig.level = c(0.0001, 0.001, 0.01, 0.05),
#  insig='label_sig',
  pch.cex=0.7, col=NULL, ncolors=200,
  wgcna_name=NULL, wgcna_name2=NULL, ...
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
  if(is.null(MEs2)){
    res <- Hmisc::rcorr(x=MEs)
  } else{

    print('here')

    # add dataset indicator to cols/rows
    d1_names <- colnames(MEs); d2_names <- colnames(MEs2);
    colnames(MEs) <- paste0(d1_names, '_D1')
    colnames(MEs2) <- paste0(d2_names, '_D2')

    res <- Hmisc::rcorr(x=MEs, y=as.matrix(MEs2))

    print('here')
    res$r <- res$r[!grepl('_D1', colnames(res$r)),grepl('_D1', colnames(res$r))]
    colnames(res$r) <- d1_names
    rownames(res$r) <- d2_names

    res$P <- res$P[!grepl('_D1', colnames(res$P)),grepl('_D1', colnames(res$P))]
    colnames(res$P) <- d1_names
    rownames(res$P) <- d2_names

  }
  res$P[is.na(res$P)] <- 0

  # plot correlogram
  corrplot::corrplot(
    res$r,
    p.mat = res$P,
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
  order_points=TRUE, restrict_range=TRUE, point_size = 0.5, alpha=1,
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

    print(cur_mod)

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

    cur_plot_df <- plot_df[,c(colnames(umap), cur_mod)]
    colnames(cur_plot_df)[3] <- "val"

    # order points:
    if(order_points == TRUE){
      cur_plot_df <- cur_plot_df %>% dplyr::arrange(val)
    } else if(order_points == "shuffle"){
      cur_plot_df <- cur_plot_df[sample(nrow(cur_plot_df)),]
    }

    # plot with ggplot
    p <- cur_plot_df %>%
      ggplot(aes_string(x=x_name, y=y_name, color="val")) +
      # ggplot(aes(x=umap1, y=umap2, color=val))
      geom_point(size=point_size, alpha=alpha) +
      ggtitle(cur_mod) + umap_theme +
      labs(color="")

    # UCell?
    if(!ucell){
      p <- p + scale_color_gradient2(
        low='grey75', mid='grey95', high=cur_color,
        breaks = plot_range,
        labels = c('-', '+'),
        guide = guide_colorbar(ticks=FALSE, barwidth=0.5, barheight=4)
      )
    } else{
      p <- p + scale_color_gradient(
        low='grey95', high=cur_color,
        breaks = plot_range,
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

  # add color to enrich_table
  mod_colors <- select(modules, c(module, color)) %>% distinct
  enrichr_df$color <- mod_colors[match(enrichr_df$module, mod_colors$module), 'color']

  # helper function to wrap text
  wrapText <- function(x, len) {
      sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
  }

  # get data to plot
  plot_df <- enrichr_df %>%
    subset(db == database & module %in% mods) %>%
    group_by(module) %>%
    top_n(n_terms, wt=Combined.Score)

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
    geom_point(aes(size=Combined.Score), color=plot_df$color) +
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
  label_center = FALSE, # only label the genes in the middle?
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

  # get hub genes:
  n_hubs <- 25
  hub_list <- lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur <- cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(n_hubs)
    colnames(cur)[2] <- 'var'
    cur %>% arrange(desc(var)) %>% .$gene_name
  })
  names(hub_list) <- mods

  print('here')

  # loop over modules
  for(cur_mod in mods){
    print(cur_mod)
    cur_color <- modules %>% subset(module == cur_mod) %>% .$color %>% unique

    # number of genes, connections
    # might make a setting to change this later but I will also have to change
    # how the graph layout works
    n_genes = 25;
    n_conns = 500;

    # name of column with current kME info
    cur_kME <- paste0('kME_', cur_mod)

    cur_genes <- hub_list[[cur_mod]]
    print(cur_genes)

    # if (length(cur_genes) < n_genes) {
    #   n_genes <-  length(cur_genes);
    #   n_conns <- n_genes * (n_genes - 1);
    # }

    # Identify the columns in the TOM that correspond to these hub genes
    matchind <- match(cur_genes, colnames(TOM))
    reducedTOM = TOM[matchind,matchind]
    orderind <- order(reducedTOM,decreasing=TRUE)

    # only  keep top connections
    connections2keep <- orderind[1:n_conns];
    reducedTOM <- matrix(0,nrow(reducedTOM),ncol(reducedTOM));
    reducedTOM[connections2keep] <- 1;

    print('here')
    print(dim(reducedTOM))
    print(n_genes)

    # only label the top 10 genes?
    if(label_center){cur_genes[11:25] <- ''}

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
      vertex.label=as.character(cur_genes),
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
    subset(!(gene_name %in% unlist(hub_list))) %>%
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

  # label vertices?

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

#' Plots the results from OverlapModulesDEGs as a bar plot
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


#' ROCCurves
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
#' ROCCurves
ROCCurves <- function(
  seurat_obj,
  roc_df=NULL,
  conf_df=NULL,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get Modules
  modules <- GetModules(seurat_obj)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get module colors:
  mod_colors <- modules %>% subset(module %in% mods) %>%
    select(c(module, color)) %>%
    distinct %>%
    arrange(module) %>% .$color

  # get the ROC info from seurat obj:
  if(is.null(roc_df) | is.null(conf_df)){
    roc_df <- GetROCData(seurat_obj, wgcna_name)$roc
    conf_df <- GetROCData(seurat_obj, wgcna_name)$conf
  }

  # plot the ROC curve
  roc_df <- roc_df %>% group_by(module) %>% arrange(sensitivity)
  conf_df <- conf_df %>% group_by(module) %>% arrange(sensitivity)
  auc_df <- distinct(roc_df[,c('module', 'auc')])

  # set factor levels for modules:
  roc_df$module <- factor(as.character(roc_df$module), levels=mods)
  conf_df$module <- factor(as.character(conf_df$module), levels=mods)
  auc_df$module <- factor(as.character(auc_df$module), levels=mods)

  p <- roc_df %>% ggplot(
    aes(x=specificity, y=sensitivity, color=module, fill=module),
  ) +
    geom_line() +
    geom_ribbon(
      data=conf_df,
      aes(x = sensitivity, ymin=lo, ymax=hi, fill=module),
      inherit.aes=FALSE, alpha=0.4
    ) +
    scale_color_manual(values = unlist(mod_colors)) +
    scale_fill_manual(values = unlist(mod_colors)) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels=c("0", "0.5", "1")) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels=c("0", "0.5", "1")) +
    xlab("1 - Specificity (FPR)") + ylab("Sensitivity (TPR)") +
    geom_text(
      data = auc_df,
      aes(color=module),
      x=0.75, y=0.1, label=paste0("AUC: ", format(auc_df$auc, digits=2)),
      inherit.aes=FALSE, size=4, color='black'
    )

  p

}



#' Displays the top n TFs in a set of modules as a bar plot
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' MotifOverlapBarPlot
MotifOverlapBarPlot <- function(
  seurat_obj,
  n_tfs = 10,
  plot_size = c(5,6),
  outdir = "MotifOverlaps/",
  motif_font = 'helvetica_regular',
  module_names = NULL,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # make output dir if it doesn't exist:
  if(!dir.exists(outdir)){dir.create(outdir)}

  # get Modules
  modules <- GetModules(seurat_obj)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  if(is.null(module_names)){module_names <- mods}

  # get overlap info from Seurat obj
  overlap_df <- GetMotifOverlap(seurat_obj, wgcna_name)
  motif_df <- GetMotifs(seurat_obj)

  # get pfm list from Seurat obj
  pfm <- GetPFMList(seurat_obj)

  # add motif ID to the overlap_df
  overlap_df$motif_ID <- motif_df$motif_ID[match(overlap_df$tf, motif_df$motif_name)]

  # subset by the modules that we are using:
  overlap_df <- overlap_df %>% subset(module %in% module_names)

  for(cur_mod in module_names){
    print(cur_mod)
    # get df for cur mod
    plot_df <- overlap_df %>% subset(module == cur_mod) %>% top_n(n_tfs, wt=odds_ratio) %>%
      arrange(desc(odds_ratio))

    # make the barplot
    p1 <- plot_df %>% ggplot(aes(y=reorder(tf, odds_ratio), fill=odds_ratio, x=odds_ratio)) +
      geom_bar(stat='identity', width=0.7) + NoLegend() +
      scale_fill_gradient(high=unique(plot_df$color), low='grey90') +
      ylab('') + theme(
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(
          t = 0, r = 0, b = 0, l = 0
        )
      )

    # make the motif logo plots:
    plot_list <- list()
    for(i in 1:nrow(plot_df)){
      cur_id <- plot_df[i,'motif_ID']
      cur_name <- plot_df[i,'tf']
      plot_list[[cur_id]] <- ggplot() +
        ggseqlogo::geom_logo( as.matrix(pfm[[cur_id]]), font=motif_font) +
        ggseqlogo::theme_logo() +
        xlab('') + ylab(cur_name) + theme(
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y = element_text(angle=0),
          plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0) # Left margin
        )
    }

    # wrap the motif logos
    patch1 <- wrap_plots(plot_list, ncol=1)

    # assemble the plot with patchwork
    outplot <- (patch1 | p1) +
          plot_layout(ncol=2, widths=c(1,2)) +
          plot_annotation(title=paste0('Motif overlaps with ', cur_mod),
          theme = theme(plot.title=element_text(hjust=0.5))
        )

    pdf(paste0(outdir, '/', cur_mod, '_motif_overlaps.pdf'), width=plot_size[1], height=plot_size[2], useDingbats=FALSE)
    print(outplot)
    dev.off()

  }

}



#' Plots gene expression of hub genes as a heatmap
#'
#' This function makes an expression heatmap of the top n hub genes per module
#' using Seurat's DoHeatmap, and then assembles them all into one big heatmap.
#'
#'
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' DoHubGeneHeatmap
DoHubGeneHeatmap <- function(
  seurat_obj,
  n_hubs = 10,
  n_cells = 200,
  group.by=NULL,
  module_names = NULL,
  combine=TRUE, #returns a list of individual heatmaps if FALSE
  draw.lines=TRUE,
  disp.min = -2.5, disp.max = 2.5, # cutoff expression values
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # use idents as grouping variable if not specified
  if(is.null(group.by)){
    group.by <- 'temp_ident'
    seurat_obj$temp_ident <- Idents(seurat_obj)
  }

  # drop if there are missing levels:
  seurat_obj@meta.data[[group.by]] <- droplevels(seurat_obj@meta.data[[group.by]])

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  modules <- modules %>% subset(module != 'grey') %>% mutate(module = droplevels(module))
  mods <- levels(modules$module)

  if(!is.null(module_names)){
    print('here')
    mods <- module_names
    modules <- modules %>% subset(module %in% mods)
  }

  # get table of module names & colors
  mod_colors <- modules %>% dplyr::select(c(module, color)) %>% distinct

  # get hub genes:
  # hub_list <- lapply(mods, function(cur_mod){
  #   cur <- subset(modules, module == cur_mod)
  #   cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
  #     top_n(n_hubs) %>% .$gene_name
  # })
  hub_list <- lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur <- cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(n_hubs)
    colnames(cur)[2] <- 'var'
    cur %>% arrange(desc(var)) %>% .$gene_name
  })
  names(hub_list) <- mods
  print(hub_list)

  seurat_obj$barcode <- colnames(seurat_obj)
  temp <- table(seurat_obj@meta.data[[group.by]])

  # sample cells
  df <- data.frame()
  for(i in 1:length(temp)){

    if(temp[[i]] < n_cells){
      cur_df <- seurat_obj@meta.data %>% subset(get(group.by) == names(temp)[i])
    } else{
      cur_df <- seurat_obj@meta.data %>% subset(get(group.by) == names(temp)[i]) %>% sample_n(n_cells);
    }
    df <- rbind(df, cur_df)
  }

  # make sampled seurat obj for plotting:
  seurat_plot <- seurat_obj %>% subset(barcode %in% df$barcode)

  plot_list <- list()
  for(i in 1:length(hub_list)){

    print(i)
    cur_mod <- names(hub_list)[i]
    print(i)
    print(hub_list[[i]])
    print(i)

    if(i == 1){
      plot_list[[i]] <- DoHeatmap(
        seurat_plot,
        features = hub_list[[i]],
        group.by=group.by,
        raster=TRUE, slot='scale.data',
        disp.min = disp.min, disp.max=disp.max,
        label=FALSE,
        group.bar=FALSE,
        draw.lines=draw.lines
      )
    } else{
      plot_list[[i]] <- DoHeatmap(
       seurat_plot,
       features=hub_list[[i]],
       group.by=group.by,
       raster=TRUE, slot='scale.data',
       group.bar.height=0,
       label=FALSE, group.bar=FALSE,
       draw.lines=draw.lines,
       disp.min = disp.min, disp.max=disp.max
     ) + NoLegend()
    }
    print(i)
    # margin:
    plot_list[[i]] <- plot_list[[i]] +
      theme(
        plot.margin = margin(0,0,0,0),
        axis.text.y = element_text(face='italic')
      )  + scale_y_discrete(position = "right")
    print(i)

  }

  # useful for ratios on colorbars
  n_total_cells <- ncol(seurat_plot)
  width_cbar <- n_total_cells / 50


  # module colorbar
  mod_colors$value <- n_hubs
  mod_colors$dummy <- 'colorbar'
  cbar_list <- list()
  for(i in 1:nrow(mod_colors)){
    cbar_list[[i]] <- mod_colors[i,] %>% ggplot(aes(y=value, x=dummy)) +
      geom_bar(position='stack', stat='identity', fill=mod_colors[i,]$color) +
      umap_theme + theme(
        plot.margin=margin(0,0,0,0)
      )
  }
  p_cbar <- wrap_plots(cbar_list, ncol=1)

  if(combine){
    out <- wrap_plots(plot_list, ncol=1) +plot_layout(guides='collect')
    out <- (p_cbar | out) + plot_layout(widths=c(width_cbar, n_total_cells))
  } else{
    out <- plot_list
  }

  out

}
