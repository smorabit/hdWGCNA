#' PlotSoftPowers
#'
#' Plot Soft Power Threshold results
#'
#' @param seurat_obj A Seurat object
#' @param selected_power power to highlight in the plots
#' @param point_size the size of the points in the plot
#' @param text_size the size of the text in the plot
#' @param plot_connectivity logical indicating whether to plot the connectivity in addition to the scale free topplogy fit.
#' @param wgcna_name The name of the WGCNA experiment in seurat_obj
#' @keywords scRNA-seq
#' @export
#' @examples
#' PlotSoftPowers
PlotSoftPowers <- function(
  seurat_obj,
  selected_power = NULL,
  point_size = 5,
  text_size=3,
  plot_connectivity = TRUE,
  wgcna_name = NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  pt <- GetPowerTable(seurat_obj, wgcna_name)

  if("group" %in% colnames(pt)){
    print("here")
    # select soft power for each group:
    power_tables <- pt %>% dplyr::group_split(group)
    soft_powers <- sapply(power_tables, function(power_table){
      power_table %>% subset(SFT.R.sq >= 0.8) %>% .$Power %>% min
    })

  } else{

    # get soft power:
    if(is.null(selected_power)){
      soft_power <- pt %>% subset(SFT.R.sq >= 0.8) %>% .$Power %>% min
    } else{
      soft_power <- selected_power
    }
    soft_powers <- NULL
    power_tables <- list("power" = pt)

  }

  plot_list <- list()

  for(i in 1:length(power_tables)){

    pt <- power_tables[[i]]

    if(!is.null(soft_powers)){
      soft_power <- soft_powers[i]
      print(i)
      print(soft_power)

      pt <- pt %>% dplyr::select(-group)
    }

    print(head(pt))

    # get other params:
    sft_r <- as.numeric(pt[pt$Power == soft_power,'SFT.R.sq'])
    mean_k <- as.numeric(pt[pt$Power == soft_power,'mean.k.'])
    median_k <- as.numeric(pt[pt$Power == soft_power,'median.k.'])
    max_k <- as.numeric(pt[pt$Power == soft_power,'max.k.'])

    # set color of text
    pt$text_color <- ifelse(
      pt$Power == soft_power, 'white', 'black'
    )

    # plot for soft power thresh:
    p1 <- pt %>% ggplot(aes(x=Power, y=SFT.R.sq)) +
      geom_rect(
        data = pt[1,],
        aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.8), fill='grey80', alpha=0.8, color=NA
      ) +
      geom_hline(yintercept = sft_r, linetype='dashed') +
      geom_vline(xintercept = soft_power, linetype= 'dashed') +
      geom_point(
        data = pt[pt$Power == soft_power,c('Power', 'SFT.R.sq')],
        aes(x=Power, y=SFT.R.sq),
        inherit.aes=FALSE,
        color = 'black',
        size=point_size
      ) +
      geom_text(label=pt$Power, color = pt$text_color, size=text_size) +
      scale_y_continuous(limits = c(0,1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      ylab('Scale-free Topology Model Fit') +
      xlab('Soft Power Threshold') +
      theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
      )

      if(plot_connectivity){

        # plot for mean connectivity:
        p2 <- pt %>% ggplot(aes(x=Power, y=mean.k.)) +
          geom_hline(yintercept = mean_k, linetype='dashed') +
          geom_vline(xintercept = soft_power, linetype= 'dashed') +
          geom_point(
            data = pt[pt$Power == soft_power,c('Power', 'mean.k.')],
            aes(x=Power, y=mean.k.),
            inherit.aes=FALSE,
            color = 'black',
            size=point_size
          ) +
          geom_text(label=pt$Power, color = pt$text_color, size=text_size) +
          scale_y_continuous(labels=scales::comma) +
          ylab('Mean Connectivity') +
          xlab('Soft Power Threshold') +
          theme(
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1)
          )

          # plot for medianan connectivity:
          p3 <- pt %>% ggplot(aes(x=Power, y=median.k.)) +
            geom_hline(yintercept = median_k, linetype='dashed') +
            geom_vline(xintercept = soft_power, linetype= 'dashed') +
            geom_point(
              data = pt[pt$Power == soft_power,c('Power', 'median.k.')],
              aes(x=Power, y=median.k.),
              inherit.aes=FALSE,
              color = 'black',
              size=point_size
            ) +
            geom_text(label=pt$Power, color = pt$text_color, size=text_size) +
            scale_y_continuous(labels=scales::comma) +
            ylab('Median Connectivity') +
            xlab('Soft Power Threshold') +
            theme(
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1)
            )

          # plot for mean connectivity:
          p4 <- pt %>% ggplot(aes(x=Power, y=max.k.)) +
            geom_hline(yintercept = max_k, linetype='dashed') +
            geom_vline(xintercept = soft_power, linetype= 'dashed') +
            geom_point(
              data = pt[pt$Power == soft_power,c('Power', 'max.k.')],
              aes(x=Power, y=max.k.),
              inherit.aes=FALSE,
              color = 'black',
              size=point_size
            ) +
            geom_text(label=pt$Power, color = pt$text_color, size=text_size) +
            scale_y_continuous(labels=scales::comma) +
            ylab('Max Connectivity') +
            xlab('Soft Power Threshold') +
            theme(
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1)
            )


          plot_list[[i]] <- list(p1,p2,p3,p4)

      } else{
        plot_list[[i]] <- p1
      }

  }

  if(length(plot_list) == 1){
    print('return')
    return(plot_list[[1]])
  }
  plot_list

}



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

    # add dataset indicator to cols/rows
    d1_names <- colnames(MEs); d2_names <- colnames(MEs2);
    colnames(MEs) <- paste0(d1_names, '_D1')
    colnames(MEs2) <- paste0(d2_names, '_D2')

    res <- Hmisc::rcorr(x=MEs, y=as.matrix(MEs2))

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
  label_legend = FALSE, ucell = FALSE, raster=FALSE, raster_dpi=500
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
      ggplot(aes_string(x=x_name, y=y_name, color="val"))

    # rasterise?
    if(raster){
      p <- p + ggrastr::rasterise(geom_point(size=point_size, alpha=alpha), dpi=raster_dpi)
    } else{
      p <- p + geom_point(size=point_size, alpha=alpha)
    }

    # add title and theme:
    p <- p + ggtitle(cur_mod) + umap_theme() + labs(color="")

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
#' Makes barplots from RunEnrhcr output.
#'
#' @param seurat_obj A Seurat object
#' @param outdir directory to place output .pdf files
#' @param n_terms the number of terms to plot in each barplot
#' @param plot_size the size of the output .pdf files (width, height)
#' @param logscale logical controlling whether to plot the enrichment on a log scale
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' EnrichrBarPlot
EnrichrBarPlot <- function(
  seurat_obj, outdir = "enrichr_plots",
  n_terms = 25, plot_size = c(6,15),
  logscale=FALSE, wgcna_name=NULL,  ...
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
#' @param database name of the enrichr database to plot.
#' @param mods names of modules to plot. All modules are plotted if mods='all' (default)
#' @param n_terms number of enriched terms to plot for each module
#' @param break_ties logical controlling whether or not to randomly select terms with equal enrichments to precisely enforce n_terms.
#' @param logscale logical controlling whether to plot the enrichment on a log scale.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' EnrichrDotPlot
EnrichrDotPlot <- function(
  seurat_obj, database, mods="all",
  n_terms = 3, break_ties=TRUE,
   logscale=TRUE, wgcna_name=NULL, ...
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
    ggplot(aes(x=module, y=Term)) +
    geom_point(aes(size=Combined.Score), color=plot_df$color) +
    RotatedAxis() +
    ylab('') + xlab('') + labs(size=lab) +
    scale_y_discrete(limits=rev) +
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
#' Visualizes the top hub genes for selected modules as a circular network plot
#'
#' @param seurat_obj A Seurat object
#' @param mods Names of the modules to plot. If mods = "all", all modules are plotted.
#' @param outdir The directory where the plots will be stored.
#' @param plot_size A vector containing the width and height of the network plots.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleNetworkPlot
ModuleNetworkPlot <- function(
  seurat_obj,
  mods="all",
  outdir="ModuleNetworks",
  plot_size = c(6,6),
  wgcna_name=NULL,
  label_center = FALSE, # only label the genes in the middle?
  edge.alpha=0.25,
  vertex.label.cex=1,
  vertex.size=6, ...
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

  # tell the user that the output is going to the output dir
  cat(paste0("Writing output files to ", outdir))

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

    # Identify the columns in the TOM that correspond to these hub genes
    matchind <- match(cur_genes, colnames(TOM))
    reducedTOM = TOM[matchind,matchind]
    orderind <- order(reducedTOM,decreasing=TRUE)

    # only  keep top connections
    connections2keep <- orderind[1:n_conns];
    reducedTOM <- matrix(0,nrow(reducedTOM),ncol(reducedTOM));
    reducedTOM[connections2keep] <- 1;

    # print('here')
    # print(dim(reducedTOM))
    # print(n_genes)

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
#' Construct a unified network plot comprising hub genes for multiple modules.
#'
#' @param seurat_obj A Seurat object
#' @param mods Names of the modules to plot. If mods = "all", all modules are plotted.
#' @param n_hubs The number of hub genes to plot for each module.
#' @param n_other The number of non-hub genes to sample from each module
#' @param edge_prop The proportion of edges in the graph to sample.
#' @param return_graph logical determining whether we return the graph (TRUE) or plot the graph (FALSE)
#' @param edge.alpha Scaling factor for the edge opacity
#' @param vertex.label.cex The font size of the gene labels
#' @param hub.vertex.size The size of the hub gene nodes
#' @param other.vertex.size The size of the other gene nodes
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' HubGeneNetworkPlot
HubGeneNetworkPlot <- function(
  seurat_obj, mods="all",
  n_hubs=6, n_other=3,
  sample_edges = TRUE,
  edge_prop = 0.5,
  return_graph=FALSE,
  edge.alpha=0.25,
  vertex.label.cex=0.5,
  hub.vertex.size=4,
  other.vertex.size=1,
  wgcna_name=NULL,
  ...
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

  print(table(selected_modules$module))

  # make sure all nodes have at least one edge!!
  edge_cutoff <- min(sapply(1:nrow(subset_TOM), function(i){max(subset_TOM[i,])}))
  edge_df <- reshape2::melt(subset_TOM) %>% subset(value >= edge_cutoff)

  edge_df$color <- future.apply::future_sapply(1:nrow(edge_df), function(i){
    gene1 = as.character(edge_df[i,'Var1'])
    gene2 = as.character(edge_df[i,'Var2'])

    col1 <- modules %>% subset(gene_name == gene1) %>% .$color
    col2 <- modules %>% subset(gene_name == gene2) %>% .$color

    if(col1 == col2){
      col = col1
    } else{
      col = 'grey90'
    }
    col
  })

  # subset edges:
  groups <- unique(edge_df$color)
  print(groups)
  if(sample_edges){
    print('here')
    # randomly sample
    temp <- do.call(rbind, lapply(groups, function(cur_group){
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_sample <- sample(1:n_edges, round(n_edges * edge_prop))
      cur_df[cur_sample,]
    }))
  } else{

    # get top strongest edges
    temp <- do.call(rbind, lapply(groups, function(cur_group){
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_df %>% dplyr::top_n(round(n_edges * edge_prop), wt=value)
    }))
  }

  edge_df <- temp
  print(dim(edge_df))

  # scale edge values between 0 and 1 for each module
  edge_df <- edge_df %>% group_by(color) %>% mutate(value=scale01(value))

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

  l <- igraph::layout_with_fr(g, ...)

  if(return_graph){return(g)}

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
#' @param return_graph logical determining whether to plot thr graph (FALSE) or return the igraph object (TRUE)
#' @param keep_grey_edges logical determining whether to show edges between genes in different modules (grey edges)
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleUMAPPlot
ModuleUMAPPlot <- function(
  seurat_obj,
  sample_edges = TRUE, # TRUE if we sample edges randomly, FALSE if we take the top edges
  edge_prop = 0.2,
  label_hubs = 5, # how many hub genes to label?
  edge.alpha=0.25,
  vertex.label.cex=0.5,
  label_genes = NULL,
  return_graph = FALSE, # this returns the igraph object instead of plotting
  keep_grey_edges = TRUE,
  wgcna_name=NULL,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get the TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # get modules,
  modules <- GetModules(seurat_obj, wgcna_name)

  # get the UMAP df:
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  mods <- levels(umap_df$modules)

  # subset module df by genes in the UMAP df:
  selected_modules <- modules[umap_df$gene,]
  selected_modules <- cbind(selected_modules, umap_df[,c('UMAP1', 'UMAP2', 'hub', 'kME')])

  # subset the TOM:
  subset_TOM <- TOM[umap_df$gene, umap_df$gene[umap_df$hub == 'hub']]

  # genes to label:
  hub_labels <- selected_modules %>% group_by(module) %>% top_n(label_hubs, wt=kME) %>% .$gene_name
  if(is.null(label_genes)){
    label_genes <- hub_labels
  } else{
    if(!any(label_genes %in% umap_df$gene)){
      stop("Some genes in label_genes not found in the UMAP.")
    }
    label_genes <- unique(c(label_genes, hub_labels))
  }
  selected_modules$label <- ifelse(selected_modules$gene_name %in% label_genes, selected_modules$gene_name, '')
  selected_modules$fontcolor <- ifelse(selected_modules$color == 'black', 'gray50', 'black')

  # set frome color
  # same color as module for all genes, black outline for the selected hub genes
  selected_modules$framecolor <- ifelse(selected_modules$gene_name %in% label_genes, 'black', selected_modules$color)

  # melt TOM into long format
  edge_df <- subset_TOM %>% reshape2::melt()
  print(dim(edge_df))

  # set color of each edge based on value:
  edge_df$color <- future.apply::future_sapply(1:nrow(edge_df), function(i){
    gene1 = as.character(edge_df[i,'Var1'])
    gene2 = as.character(edge_df[i,'Var2'])

    col1 <- selected_modules[selected_modules$gene_name == gene1, 'color']
    col2 <- selected_modules[selected_modules$gene_name == gene2, 'color']

    if(col1 == col2){
      col = col1
    } else{
      col = 'grey90'
    }
    col
  })

  # keep grey edges?
  if(!keep_grey_edges){
    edge_df <- edge_df %>% subset(color != 'grey90')
  }

  # subset edges:
  groups <- unique(edge_df$color)
  if(sample_edges){
    # randomly sample
    temp <- do.call(rbind, lapply(groups, function(cur_group){
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_sample <- sample(1:n_edges, round(n_edges * edge_prop))
      cur_df[cur_sample,]
    }))
  } else{

    # get top strongest edges
    temp <- do.call(rbind, lapply(groups, function(cur_group){
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_df %>% dplyr::top_n(round(n_edges * edge_prop), wt=value)
    }))
  }

  edge_df <- temp
  print(dim(edge_df))

  # scale edge values between 0 and 1 for each module
  edge_df <- edge_df %>% group_by(color) %>% mutate(value=scale01(value))

  # edges & vertices are plotted in igraph starting with the first row, so re-order s.t. strong edges are on bottom, all gray on the top of the table:
  edge_df <- edge_df %>% arrange(value)
  edge_df <- rbind(
    subset(edge_df, color == 'grey90'),
    subset(edge_df, color != 'grey90')
  )

  # set alpha of edges based on kME
  edge_df$color_alpha <- ifelse(
    edge_df$color == 'grey90',
    alpha(edge_df$color, alpha=edge_df$value/2),
    alpha(edge_df$color, alpha=edge_df$value)
  )

  # re-order vertices so hubs are plotted on top
  selected_modules <- rbind(
    subset(selected_modules , hub == 'other'),
    subset(selected_modules , hub != 'other')
  )

  # re-order vertices so labeled genes are on top
  selected_modules <- rbind(
    subset(selected_modules , label == ''),
    subset(selected_modules , label != '')
  )

  # setup igraph:
  g <- igraph::graph_from_data_frame(
    edge_df,
    directed=FALSE,
    vertices=selected_modules
  )

  if(return_graph){return(g)}

  plot(
    g,
    layout=  as.matrix(selected_modules[,c('UMAP1', 'UMAP2')]),
    # edge.color=adjustcolor(E(g)$color, alpha.f=edge.alpha),
    edge.color=adjustcolor(E(g)$color_alpha, alpha.f=edge.alpha),
    vertex.size=V(g)$kME * 3,
    edge.curved=0,
    edge.width=0.5,
    vertex.color=V(g)$color,
    vertex.label=V(g)$label,
    vertex.label.dist=1.1,
    vertex.label.degree=-pi/4,
    vertex.label.family='Helvetica', #vertex.label.font=vertex_df$font,
    vertex.label.font = 3,
    vertex.label.color = V(g)$fontcolor,
    vertex.label.cex=0,
    vertex.frame.color=V(g)$framecolor,
    margin=0
  )

}


#' OverlapDotPlot
#'
#' Makes barplots from Enrichr data
#'
#' @param overlap_df the Module/DEG overlap table from OverlapModulesDEGs
#' @param plot_var the name of the overlap statistic to plot
#' @param logscale logical controlling whether to plot the result on a log scale, useful for odds ratio
#' @param neglog logical controlling wehether to plot the result as a negative log, useful for p-value / FDR
#' @param plot_significance logical controlling whether to plot the significance levels on top of the dots
#' @keywords scRNA-seq
#' @export
#' @examples
#' OverlapDotPlot
OverlapDotPlot <- function(
  overlap_df, plot_var = 'odds_ratio',
  logscale=TRUE,
  neglog=FALSE,
  plot_significance=TRUE,
  ...
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
#' Plots the results from OverlapModulesDEGs as a bar plot
#'
#' @param overlap_df the Module/DEG overlap table from OverlapModulesDEGs
#' @param plot_var the name of the overlap statistic to plot
#' @param logscale logical controlling whether to plot the result on a log scale, useful for odds ratio
#' @param neglog logical controlling wehether to plot the result as a negative log, useful for p-value / FDR
#' @param label_size the size of the module labels in the bar plot
#' @keywords scRNA-seq
#' @export
#' @examples
#' OverlapBarPlot
OverlapBarPlot <- function(
  overlap_df,
  plot_var = 'odds_ratio',
  logscale=FALSE, neglog=FALSE,
  label_size=2,
  ...
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
    overlap_df[[plot_var]] <- -1 * log(overlap_df[[plot_var]])
    label <- paste0('-log(', label, ')')
    yint = -1 * log(yint)
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

      # add the labels:
      p <- p +
        geom_text(
          aes(label=module, x=module, y=get(plot_var)), color='black', size=label_size, hjust='inward'
        )

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
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
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
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
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
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
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
      umap_theme() + theme(
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

#' PlotModulePreservatgion
#'
#' Plotting function for Module Preservation statistics
#'
#' @param seurat_obj A Seurat object
#' @param name The name of the module preservation analysis to plot given in ModulePreservation
#' @param statistics Which module preservation statistics to plot? Choices are summary, all, or a custom list
#' @param plot_labels logical determining whether to plot the module labels
#' @param label_size the size of the module labels
#' @param mod_point_size the size of the points in each plot
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' PlotModulePreservation
PlotModulePreservation <- function(
  seurat_obj,
  name,
  statistics = 'summary', # can be summary, rank, all, or a custom list
  plot_labels = TRUE,
  label_size = 4,
  mod_point_size = 4,
  wgcna_name = NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get the module preservation stats:
  mod_pres <- GetModulePreservation(seurat_obj, name, wgcna_name)
  obs_df <- mod_pres$obs
  Z_df <- mod_pres$Z

  # get module colors:
  modules <- GetModules(seurat_obj, wgcna_name)
  module_colors <- modules %>% dplyr::select(c(module, color)) %>% distinct
  mods <- rownames(Z_df)
  mod_colors <- module_colors$color[match(mods, module_colors$module)]
  mod_colors = ifelse(is.na(mod_colors), 'gold', mod_colors)

  # what are we going to plot?
  if(statistics == 'summary'){
    stat_list <- c("Zsummary.qual", "Zsummary.pres")
  } else if(statistics == 'rank'){
    stat_list <- colnames(obs_df[,-1])[grepl("Rank", colnames(obs_df[,-1]))]
  } else if(statistics == 'all'){
    stat_list <- c(colnames(obs_df[,-1])[grepl("Rank", colnames(obs_df[,-1]))], colnames(Z_df[,-1]))
  } else{
    stat_list <- statistics
  }

  stat_list <- stat_list[stat_list != 'moduleSize']

  plot_list <- list()
  for(statistic in stat_list){

    print(statistic)

    if(statistic %in% colnames(obs_df)){
      values <- obs_df[,statistic]
    } else if(statistic %in% colnames(Z_df)){
      values <- Z_df[,statistic]
    } else{
      stop("Invalid name for statistic.")
    }

    # setup plotting df
    plot_df <- data.frame(
      module = mods,
      color = mod_colors,
      value = values,
      size = Z_df$moduleSize
    )

    # don't include grey & gold:
    plot_df <- plot_df %>% subset(!(module %in% c('grey', 'gold')))


    if(grepl("Rank", statistic)){
      cur_p <-  plot_df %>% ggplot(aes(x=size, y=value, fill=module, color=module)) +
        geom_point(size=mod_point_size, pch=21, color='black') +
        scale_y_reverse()
    } else{
      cur_p <- plot_df %>% ggplot(aes(x=size, y=value, fill=module, color=module)) +
        geom_rect(
          data = plot_df[1,],
          aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=2), fill='grey75', alpha=0.8, color=NA) +
        geom_rect(
          data=plot_df[1,],
          aes(xmin=0, xmax=Inf, ymin=2, ymax=10), fill='grey92', alpha=0.8, color=NA) +
        geom_point(size=mod_point_size, pch=21, color='black')
    }

    cur_p <- cur_p +
      scale_fill_manual(values=plot_df$color) +
      scale_color_manual(values=plot_df$color) +
      scale_x_continuous(trans='log10') +
      ylab(statistic) +
      xlab("Module Size") +
      ggtitle(statistic) +
      NoLegend() +
      theme(
        plot.title = element_text(hjust = 0.5)
      )


    if(plot_labels){
      cur_p <- cur_p + geom_text_repel(label = plot_df$module, size=label_size)
    }

    plot_list[[statistic]] <- cur_p

  }

  if(length(plot_list) == 1){return(plot_list[[1]])}

  plot_list

}

#' PlotModuleTraitCorrelation
#'
#' Plotting function for Module Preservation statistics
#'
#' @param seurat_obj A Seurat object
#' @param
#' @param
#' @param
#' @param
#' @param plot_labels logical determining whether to plot the module labels#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' PlotModulePreservation
PlotModuleTraitCorrelation <- function(
  seurat_obj,
  high_color = 'red',
  mid_color = 'grey90',
  low_color = 'blue',
  label = NULL,
  label_symbol = 'stars',
  plot_max = NULL,
  text_size = 2,
  text_color = 'black',
  text_digits = 3,
  combine = TRUE,
  wgcna_name = NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}


  # get the module trait correlation results:
  temp <- GetModuleTraitCorrelation(seurat_obj)
  cor_list <- temp$cor
  pval_list <- temp$pval
  fdr_list <- temp$fdr

  # get module colors:
  modules <- GetModules(seurat_obj, wgcna_name)
  module_colors <- modules %>%
    dplyr::select(c(module, color)) %>%
    distinct %>% subset(module != 'grey') %>%
    arrange(module)
  mod_colors <- module_colors$color

  # dummy variable
  module_colors$var <- 1

  # make the colorbar as its own heatmap
  module_colorbar <- module_colors %>%
    ggplot(aes(x=module, y=var, fill=module)) +
    geom_tile() +
    scale_fill_manual(values=mod_colors) +
    NoLegend() +
    RotatedAxis() +
    theme(
      plot.title=element_blank(),
      axis.line=element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      plot.margin=margin(0,0,0,0)
    )


  plot_list <- list()
  for(i in names(cor_list)){
    cor_mat <- as.matrix(cor_list[[i]])
    pval_mat <- as.matrix(pval_list[[i]])
    fdr_mat <- as.matrix(fdr_list[[i]])
    print(i)

    plot_df <- reshape2::melt(cor_mat)
    colnames(plot_df) <- c("Trait", "Module", "cor")

    #p_df <- reshape2::melt(pval_mat)
    if(!is.null(label)){
      if(label == 'fdr'){
        p_df <- reshape2::melt(fdr_mat)
      } else if(label == 'pval'){
        p_df <- reshape2::melt(pval_mat)
      }
      colnames(p_df) <- c("Trait", "Module", "pval")

      # add pval to plot_df
      plot_df$pval <- p_df$pval
      print(levels(plot_df$Trait))

      if(label_symbol == 'stars'){
        plot_df$significance <- gtools::stars.pval(plot_df$pval)
      } else if(label_symbol == 'numeric'){
        plot_df$significance <- ifelse(
          plot_df$pval <= 0.05,
          formatC(plot_df$pval, digits=text_digits), ''
        )
      } else{
        stop('Invalid input for label_symbol. Valid choices are stars or numeric.')
      }
    }

    # get limits for plot:
    if(is.null(plot_max)){
      max_plot <- max(abs(range(plot_df$cor)))
    } else{
      max_plot <- plot_max

      # fix values outside of the specified range:
      plot_df$cor <- ifelse(abs(plot_df$cor) >= plot_max, plot_max * sign(plot_df$cor), plot_df$cor)
    }

    p <- ggplot(plot_df, aes(x=Module, y=as.numeric(Trait), fill=cor)) +
      geom_tile() +
      scale_fill_gradient2(
        limits=c(-1*max_plot,max_plot),
        high=high_color,
        mid=mid_color,
        low=low_color,
        guide = guide_colorbar(ticks=FALSE, barwidth=16, barheight=0.5)
      ) +
      scale_y_continuous(
        breaks = 1:length(levels(plot_df$Trait)),
        labels=levels(plot_df$Trait),
        sec.axis = sec_axis(
          ~.,
          breaks = 1:length(levels(plot_df$Trait)),
          labels=levels(plot_df$Trait)
        )
      ) +
      RotatedAxis() + ylab('') + xlab('') + ggtitle(i) +
      # labs(fill = 'Correlation') +
      theme(
        plot.title=element_text(hjust=0.5),
        axis.line=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y.left = element_blank(),
        legend.title = element_blank(),
        legend.position='bottom'
      )

    if(!is.null(label)){
      p <- p + geom_text(label=plot_df$significance, color=text_color, size=text_size)
    }

    plot_list[[i]] <- p

  }

  if(combine){

    #plot_list <- c(plot_list, cbar)
    #names(plot_list)[length(plot_list)] <- 'module'

    for(i in 1:length(plot_list)){

      plot_list[[i]] <- plot_list[[i]] +
        ylab(names(plot_list)[i]) +
        theme(
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          legend.position='bottom',
          axis.text.x = element_blank(),
          axis.ticks=element_blank()
        )
    }

    # assemble with patchwork:
    plot_list[['module']] <- module_colorbar

    out <- wrap_plots(plot_list, ncol=1) +
      plot_layout(
        guides = 'collect',
        heights = c(rep(1, length(plot_list)-1), 0.15)
      ) +
      plot_annotation(
        theme=theme(
          plot.title=element_text(hjust=0.5),
          legend.position = 'bottom',
          legend.justification = 0.5
        )
      )

    return(out)

  } else{
    return(plot_list)
  }

}



#' ModuleTFNetwork
#'
#' Plotting the relationships between a TF and the co-expression modules
#'
#' @param seurat_obj A Seurat object
#' @param tf_name the Motif name for the tf
#' @param tf_gene_name the gene associated with this tf in the rownames(seurat_obj)
#' @param edge.alpha scaling factor for edge opacity in the network
#' @param cor_thresh threshold to plot correlation edges between modules
#' @param high_color color for positive correlation
#' @param mid_color color for zero correlation
#' @param low_color color for negative correlation
#' @param slot the slot in the seurat object to extract expression data for the tf_gene_name
#' @param size.scale scaling factor for the size of each node
#' @param tf_x x coordinate for the TF if the TF is not found in the UMAP
#' @param tf_y y coordinate for the TF if the TF is not foudn in the UMAP
#' @param wgcna_name the name of the WGCNA experiment in the seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleTFNetwork
ModuleTFNetwork <- function(
  seurat_obj,
  tf_name,
  tf_gene_name,
  edge.alpha = 0.75,
  cor_thresh =  0.25,
  high_color = 'red',
  mid_color = 'grey',
  low_color = 'blue',
  slot = 'data',
  size.scale = 30,
  tf_x = 0,
  tf_y = 0,
  wgcna_name = NULL

){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules,
  modules <- GetModules(seurat_obj, wgcna_name) %>%
    subset(module != 'grey') %>%
    mutate(module = droplevels(module))
  MEs <- GetMEs(seurat_obj, TRUE, wgcna_name) %>% as.matrix
  MEs <- MEs[,colnames(MEs) != 'grey']
  mod_sizes <- table(modules$module)

  # correlation of MEs:
  module_cor <- Hmisc::rcorr(x=MEs, type='pearson')$r
  module_cor[lower.tri(module_cor)] <- NA
  module_cor <- reshape2::melt(module_cor) %>% na.omit
  module_cor <- subset(module_cor, abs(value) >= cor_thresh & Var1 != Var2)
  module_cor

  # correlation of modules with TF expression:
  cur_exp <- GetAssayData(seurat_obj, slot=slot)[tf_gene_name,]
  exp_cor <- Hmisc::rcorr(x=MEs, y=cur_exp)$r[1:ncol(MEs),'y']
  exp_cor <- data.frame(
    mod = names(exp_cor),
    value = as.numeric(exp_cor)
  )

  plot_lim <- abs(max(c(abs(range(exp_cor$value)), abs(range(module_cor$value)))))

  # make a dummy ggplot so I can get the colors:
  p <- ggplot(module_cor, aes(x=Var1, y=Var2, color=value)) +
    geom_point() +
    scale_color_gradient2(high=high_color, mid=mid_color, low=low_color, limits=c(-1*plot_lim, plot_lim))
  ggp <- ggplot_build(p)
  module_cor$color <- ggp$data[[1]]$colour

  p <- ggplot(exp_cor, aes(x=mod, y=mod, color=value)) +
    geom_point() +
    scale_color_gradient2(high=high_color, mid=mid_color, low=low_color, limits=c(-1*plot_lim, plot_lim))
  ggp <- ggplot_build(p)
  exp_cor$color <- ggp$data[[1]]$colour

  # get the UMAP df:
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  mods <- levels(umap_df$modules)

  # get tf info
  tf_match <- GetMotifMatrix(seurat_obj)
  tf_targets <- GetMotifTargets(seurat_obj)
  motif_df <- GetMotifs(seurat_obj)
  overlap_df <- GetMotifOverlap(seurat_obj)

  # compute the UMAP centroids:
  centroid_df <-
    umap_df %>% dplyr::select(c(UMAP1, UMAP2, module)) %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

  # subset overlap_df by current TF:
  cur_overlap <- subset(overlap_df, tf == tf_name)

  # combine the two datasets:
  node_df <- dplyr::left_join(centroid_df, cur_overlap, by='module') %>%
    dplyr::rename(c(UMAP1=x, UMAP2=y, name=module))

  node_df$size <- as.numeric(node_df$size_intersection) / as.numeric(mod_sizes)

  # add info for tf to this table:
  if(tf_gene_name %in% umap_df$gene){
    tf_df <- umap_df[umap_df$gene == tf_gene_name, ] %>%
      dplyr::select(-c(hub, kME, module)) %>%
      dplyr::rename(name=gene)
  } else{
    tf_df <- data.frame(
      name = tf_gene_name,
      module = 'grey',
      color = 'grey',
      UMAP1 = tf_x,
      UMAP2 = tf_y
    )
  }

  tf_df$size <- 0.25

  # set up edge df
  edge_df <- data.frame(
    Var1 = tf_gene_name,
    Var2 = node_df$name,
    value = node_df$odds_ratio,
    color = exp_cor$color
  )

  # set the edge color to grey if the overlap isn't significant:
  edge_df$color <- ifelse(cur_overlap$fdr <= 0.05, edge_df$color, 'grey')

  # set up node df
  node_df <- dplyr::bind_rows(node_df, tf_df) %>% as.data.frame()
  rownames(node_df) <- as.character(node_df$name)

  # set up directed edge df:
  g1 <- igraph::graph_from_data_frame(
    edge_df,
    directed=TRUE,
    vertices=node_df
  )

  # set up undirected net
  g2 <- igraph::graph_from_data_frame(
    module_cor,
    directed=FALSE,
    vertices=node_df
  )


  plot(
    g2,
    layout = as.matrix(node_df[,c('UMAP1', 'UMAP2')]),
    vertex.size=1,
    edge.curved=0,
    edge.width=1,
    vertex.color='grey',
    vertex.label='',
    edge.color=adjustcolor(E(g2)$color, alpha.f=edge.alpha),
  )

  plot(
    g1,
    layout = as.matrix(node_df[,c('UMAP1', 'UMAP2')]),
    edge.color=adjustcolor(E(g1)$color, alpha.f=edge.alpha),
    vertex.size=V(g1)$size * size.scale,
    edge.curved=0,
    edge.width=edge_df$value*2,
    vertex.color=V(g1)$color,
    vertex.label=V(g1)$name,
    vertex.label.dist=1.1,
    vertex.label.degree=-pi/4,
    vertex.label.family='Helvetica',
    vertex.label.font = 3,
    vertex.label.color = 'black',
    vertex.label.cex=0,
    vertex.frame.color='black',
    margin=0,
    edge.arrow.size=edge_df$value/2,
    add=TRUE
  )


}

PlotKMEs <- function(
  seurat_obj,
  n_hubs=10,
  text_size=2,
  ncol = 5,
  plot_widths = c(3,2),
  wgcna_name = NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 'grey')
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  mod_colors <- modules %>% subset(module %in% mods) %>%
    select(c(module, color)) %>%
    distinct


  #get hub genes:
  hub_df <- do.call(rbind, lapply(mods, function(cur_mod){
    print(cur_mod)
    cur <- subset(modules, module == cur_mod)
    cur <- cur[,c('gene_name', 'module', paste0('kME_', cur_mod))]
    names(cur)[3] <- 'kME'
    cur <- dplyr::arrange(cur, kME)
    top_genes <- cur %>% dplyr::top_n(n_hubs, wt=kME) %>% .$gene_name
    cur$lab <- ifelse(cur$gene_name %in% top_genes, cur$gene_name, "")
    cur
  }))
  head(hub_df)

  plot_list <- lapply(mods, function(x){
    print(x)
    cur_color <- subset(mod_colors, module == x) %>% .$color
    cur_df <- subset(hub_df, module == x)
    top_genes <- cur_df %>% dplyr::top_n(n_hubs, wt=kME) %>% .$gene_name
    p <- cur_df %>% ggplot(aes(x = reorder(gene_name, kME), y = kME)) +
      geom_bar(stat='identity', width=1, color = cur_color, fill=cur_color) +
      ggtitle(x) +
      #xlab(paste0('kME_', x)) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.title.x = element_blank(),
        axis.line.x = element_blank()
      )
    p_anno <- ggplot() + annotate(
      "label",
      x = 0,
      y = 0,
      label = paste0(top_genes, collapse="\n"),
      size=text_size,
      fontface = 'italic',
      label.size=0
    ) + theme_void()
    patch <- p + p_anno + plot_layout(widths=plot_widths)
    patch
  })

  wrap_plots(plot_list, ncol=ncol)

}
