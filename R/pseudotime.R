
#' BinPseudotime
#'
#' Makes evenly-spaced bins of cells along a pseudotime trajectory.
#'
#' @param seurat_obj A Seurat object
#' @param n_bins number of hub genes to use in the UMAP computation
#' @param pseudotime_col name of the column in seurat_obj@meta.data containing pseudotime information
#' @export
BinPseudotime <- function(
  seurat_obj,
  n_bins = 50,
  pseudotime_col = 'pseudotime'
){

  # cut pseudotime into bins
  for(cur_bins in n_bins){
    bin_name <- paste0(pseudotime_col, '_bins_', cur_bins)
    seurat_obj@meta.data[[bin_name]] <- as.numeric(seurat_obj@meta.data[[pseudotime_col]]) %>% Hmisc::cut2(g=cur_bins)

    # make slot in the seurat misc related to these bins:
    seurat_obj@misc[[bin_name]] <- list()

  }

  # return seurat obj
  seurat_obj
}

#' PlotModuleTrajectory
#'
#' Plots module eigengene dynamics over a pseudotime trajectory.
#'
#' @param seurat_obj A Seurat object
#' @param pseudotime_col name of the column in seurat_obj@meta.data containing pseudotime information. Multiple names can be passed to plot multiple trajectories simultaneously.
#' @param n_bins number of pseudotime bins/windows used to smooth the MEs
#' @param harmonized logical determining whether or not to use the harmonized MEs
#' @param point_size size of the points in the plot
#' @param line_size width of the lines in the plot
#' @param n_col number of columns for different modules in the plot
#' @param se logical determining whether or not to show the standard error on the plot
#' @param group_colors optional list of colors to differentiate multiple pseudotime trajectories. Must be the same length as pseudotime_col
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @import ggplot2
#' @import Seurat
#' @import patchwork
#' @export
PlotModuleTrajectory <- function(
    seurat_obj,
    pseudotime_col = 'pseudotime',
    n_bins = 20,
    harmonized=TRUE,
    ncol = 4, 
    point_size = 1,
    line_size = 1,
    se = TRUE,
    group_colors = NULL,
    wgcna_name = NULL
){

    # Note: if length(pseudotime_col) > 1, make the multi-trajectory plot
    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

    # check for the group_colors:
    if(!is.null(group_colors)){
        if(length(group_colors) != length(pseudotime_col)){
            stop('group_colors should contain the same number of elements as pseudotime_col')
        }
        group_cp <- group_colors; names(group_cp) <- pseudotime_col
    }

    # run the pseudotime binning:
    for(ps_col in pseudotime_col){
        seurat_obj <- BinPseudotime(
            seurat_obj, 
            pseudotime_col=ps_col, 
            n_bins=n_bins
        )
    }

    # get relevant module information
    MEs <- GetMEs(seurat_obj, harmonized=harmonized, wgcna_name=wgcna_name)
    modules <- GetModules(seurat_obj, wgcna_name=wgcna_name)
    mods <- levels(modules$module)
    mods <- mods[mods!='grey']
    module_colors <- modules %>% dplyr::select(c(module, color)) %>% dplyr::distinct()
    
    #rownames(module_colors) <- module_colors$module
    mod_colors <- module_colors$color; names(mod_colors) <- module_colors$module

    # merge the MEs with the seurat metadata
    plot_df <- cbind(seurat_obj@meta.data, MEs)

    print(head(plot_df))

    # loop through each pseudotime trajectory:
    avg_list <- list()
    for(ps_col in pseudotime_col){

        ps_bin_name <- paste0(ps_col, '_bins_', n_bins)

        # compute the average ME in each pseudotime bin
        avg_scores <- plot_df %>%
                group_by(get(ps_bin_name)) %>%
                select(all_of(mods)) %>%
                summarise_all(mean)

        colnames(avg_scores)[1] <- 'bin'
        avg_df <- reshape2::melt(avg_scores)
        avg_df$bin <- as.numeric(avg_df$bin)
        avg_df$group <- ps_col
        avg_list[[ps_col]] <- avg_df
    }


    plot_df <- do.call(rbind, avg_list)

    # plotting a single trajectory (color by module)
    if(length(pseudotime_col) == 1){
        p <- ggplot(
            plot_df,
            aes(x = bin, y=value, color=variable, fill=variable)
        ) + 
        geom_point(size=point_size) + 
        geom_smooth(size=line_size, se=se) + 
        geom_hline(yintercept=0, linetype='dashed', color='grey') +
        scale_color_manual(values=mod_colors) +
        scale_fill_manual(values=mod_colors) +
        xlab('Pseudotime') + 
        ylab('Module Eigengene') + 
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(size=1, fill=NA, color='black')
        ) +
        NoLegend()


    } else{

        # plotting multiple trajectories (color by group)

        p <- ggplot(
            plot_df,
            aes(x = as.numeric(bin), y = value, color=group)
        ) +
        geom_smooth(se=FALSE, size=line_size) +
        geom_hline(yintercept=0, linetype='dashed', color='grey') +
        xlab('Pseudotime') +
        ylab('Module Eigengene') +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
               axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(size=1, fill=NA, color='black')
        ) +
        labs(color = 'Trajectory')

        if(!is.null(group_colors)){
            p <- p + scale_color_manual(values=group_cp)

        }

    }

    # assemble with patchwork
    patch <- p + facet_wrap(~variable, ncol=ncol, scales='free')
    patch

}
