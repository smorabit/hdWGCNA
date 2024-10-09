#' RegulonBarPlot
#'
#' Plots the top target genes within a specific TF regulon
#'
#' @return ggplot object containing the RegulonBarPlot
#'
#' @param seurat_obj A Seurat object
#' @param selected_tfs A TF whose target genes will be shown
#' @param cutoff Cutoff for the edge weights, links below this value will not be included.
#' @param top_n Number of top and bottom target genes to include in the plot.
#' @param TFs_only Logical indicating whether or not to use only TFs in the plot or to use TFs and other genes.
#' @param high_color Color corresponding to positive edge weights
#' @param mid_color Color corresponding to edge weights close to 0
#' @param low_color Color corresponding to negative edge weights
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' RegulonBarPlot creates a bar plot showing the top target genes for a particular TF 
#' based on the TF regulatory network analysis. Target genes are ranked by the predicted 
#' interaction strength from XGBoost (Gain) multiplied by the sign of the correlation 
#' between the TF and the target gene. 
#' 
#' @import Seurat
#' @import ggplot2
#' @export
RegulonBarPlot <- function(
    seurat_obj, 
    selected_tf, 
    cutoff=0.2,
    top_n=Inf,
    TFs_only = FALSE,
    high_color = 'orange2',
    mid_color = 'white',
    low_color = 'dodgerblue',
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # get the regulon df
    tf_regulons <- GetTFRegulons(seurat_obj, wgcna_name)
    if(TFs_only){tf_regulons <- subset(tf_regulons, gene %in% unique(tf_regulons$tf))}

    if(! selected_tf %in% tf_regulons$tf){
        stop(paste0('Invalid selection for selected_tf. Make sure you select a TF that is present in the dataset.'))
    }

    # set up dataframe for plotting
    plot_df <- tf_regulons %>% 
        subset(tf == selected_tf) %>% 
        mutate(score = Gain*sign(Cor)) %>%
        subset(abs(score) > cutoff) 

    top <- plot_df %>% subset(Cor > 0) %>% slice_max(n=top_n, order_by=Gain)
    bottom <- plot_df %>% subset(Cor < 0) %>% slice_min(n=top_n, order_by=Gain)
    plot_df <- rbind(top, bottom)

    # make the barplot
    p <- plot_df %>%
        ggplot(aes(x=score, y=reorder(gene, score), fill=score)) + 
        geom_bar(stat='identity', width=1) + 
        scale_fill_gradient2(low=low_color, mid=mid_color, high=high_color) + 
        geom_vline(xintercept=0, color='black') +
        geom_text(aes(label=reorder(gene, score)), color='black', size=3.5,  fontface='italic', hjust='inward') +
        ggtitle(paste0(selected_tf, ' predicted targets')) +
        xlab('Regulatory score') +
        theme(
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            plot.title = element_text(hjust=0.5),
            axis.title.y = element_blank(),
            axis.line.y = element_blank(),
        ) + Seurat::NoLegend()

    p

}