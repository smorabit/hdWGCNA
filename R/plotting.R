
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

#' MECorrelogram
#'
#' Plot Module Eigengene correlogram
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' MECorrelogram
MECorrelogram <- function(
  seurat_obj, wgcna_name = NULL, exclude_grey=TRUE,
  type = 'upper', order='original', method='ellipse',
  tl.col = 'black', tl.srt=45,
  sig.level = c(0.0001, 0.001, 0.01, 0.05),
#  insig='label_sig',
  pch.cex=0.7, col=NULL, ncolors=200, ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get MEs and convert to matrix
  MEs <- as.matrix(GetMEs(seurat_obj, wgcna_name))

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
  reduction='umap', cor_cutoff=0.2, label_vertices=FALSE, edge_scale=5,
  vertex_size=15, niter=100, vertex_frame=FALSE
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # Get module eigengenes
  MEs <- GetMEs(seurat_obj, wgcna_name)
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

  # plot colorbar:
  # TODO
  # Possibly get rid of this because it depends on a super random package
  # (fields) and kinda looks ugly?
  # colfunc <- colorRampPalette(c( "seagreen", "white", "darkorchid1"))
  # image.plot(legend.only=T, zlim=c(-1,1), col=colfunc(256) )
  #

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





#' PlotEnrichr
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
#' PlotEnrichr
PlotEnrichr <- function(
  seurat_obj, outdir = "enrichr_plots",
  n_terms = 25, plot_size = c(6,15),
  wgcna_name=NULL, logscale=FALSE, ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- as.character(unique(modules$module))
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
        plot_df$Combined.Score <- log2(plot_df$Combined.Score)
        lab <- 'Enrichment log2(combined score)'
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
