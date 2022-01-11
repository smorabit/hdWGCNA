
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
  insig='label_sig', pch.cex=0.7, col=NULL, ncolors=200, ...){

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
    insig=insig, pch.cex=pch.cex,
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
  if(is.null(clusters)){
    cluster_col <- Idents(seurat_obj)
  } else{
  }
  clusters <- droplevels(seurat_obj$annotation)
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
  l <- igraph::layout_with_fr(
    g, weights=E(g)$value,
    coords = as.matrix(data.frame(x=V(g)$x, y=V(g)$y)),
    niter=niter
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
  colfunc <- colorRampPalette(c( "seagreen", "white", "darkorchid1"))
  image.plot(legend.only=T, zlim=c(-1,1), col=colfunc(256) )


}


#' MEFeaturePlot
#'
#' Plot module eigengenes as a FeaturePlot
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' MEFeaturePlot
MEFeaturePlot<- function(
  seurat_obj, module_names=NULL, wgcna_name = NULL,
  harmonized=TRUE, reduction='umap',
  order=TRUE, restrict_range=TRUE, point_size = 0.5,
  label_legend = FALSE
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get MEs, module data from seurat object
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)

  # use all modules except gray if not specified by the user
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
    if(order){
      plot_df <- plot_df %>% dplyr::arrange_(cur_mod)
    }

    # label for plot:
    if(harmonized){
      label <- paste0('hME_', cur_mod)
    } else{
      label <- paste0('ME_', cur_mod)
    }

    # plot with ggplot
    plot_list[[cur_mod]] <- plot_df %>%
      ggplot(aes_string(x=x_name, y=y_name, color=cur_mod)) +
      geom_point(size=point_size) +
      scale_color_gradient2(
        low='grey75', mid='grey95', high=cur_color,
        breaks = c(plot_range[1], 0, plot_range[2]),
        labels = c('-', '0', '+'),
      ) +
      ggtitle(label) + umap_theme +
      labs(color="")
  }

  # return plot
  if(length(plot_list) == 1){
    p <- plot_list[[1]]
  } else{
    p <- plot_list
  }

  p

}
