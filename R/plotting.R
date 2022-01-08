
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
  seurat_obj, groupLabels="Module colors",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = "", ...
){

  WGCNA::plotDendroAndColors(
    seurat_obj@misc$wgcna_net$dendrograms[[1]],
    as.character(seurat_obj@misc$wgcna_net$colors),
    groupLabels=groupLabels,
    dendroLabels = dendroLabels,
    hang = hang,
    addGuide = addGuide,
    guideHang = guideHang,
    main = main,
    ...
  )
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
  seurat_obj, modules,
  harmonized=TRUE, reduction='umap',
  order=TRUE, restrict_range=TRUE, point_size = 0.5,
  label_legend = FALSE
){

  # get MEs from seurat object
  MEs <- GetMEs(seurat_obj, harmonized)

  # get  reduction from seurat obj
  umap <- seurat_obj@reductions[[reduction]]@cell.embeddings
  x_name <- colnames(umap)[1]
  y_name <- colnames(umap)[2]

  # merge into one df for plotting
  plot_df <- cbind(umap, MEs) %>% as.data.frame()

  plot_list <- list()
  for(cur_mod in modules){

    cur_color <- cur_mod

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
      geom_point(size=0.5) +
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
