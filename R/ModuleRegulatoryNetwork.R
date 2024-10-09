


#' ModuleRegulatoryNetwork
#'
#' Summarizes Transcription Factor regulatory networks across co-expression modules.
#'
#' @param seurat_obj A Seurat object containing the single-cell data and WGCNA results.
#' @param TFs_only Logical; if TRUE (default), only transcription factor (TF) genes are 
#' included in the regulatory network. If FALSE, the network includes all genes.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#'
#' @details
#' This function summarizes transcription factor regulatory networks across modules to infer module-module relationships. The number of directed TF links across modules are counted and normalized. The strength of these links from the XGBoost model are also tracked.
#'
#' @return A data frame with the following columns:
#' \item{source}{The source module in the regulatory interaction.}
#' \item{target}{The target module in the regulatory interaction.}
#' \item{n_pos}{The number of positive regulatory links from the source module to the target module.}
#' \item{sum_pos}{The sum of the gain values for the positive links.}
#' \item{n_neg}{The number of negative regulatory links from the source module to the target module.}
#' \item{sum_neg}{The sum of the gain values for the negative links.}
#' \item{mean_pos}{The average gain of the positive links.}
#' \item{mean_neg}{The average gain of the negative links.}
#' \item{score_pos}{The number of positive links normalized by the number of TFs in the source module.}
#' \item{score_neg}{The number of negative links normalized by the number of TFs in the source module.}
#'
#' @import Seurat 
#' @import igraph 
#' @import ggraph 
#' @import tidygraph
#' @export
ModuleRegulatoryNetwork <- function(
    seurat_obj,
    TFs_only = TRUE,
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # get regulon & TF net info
    tf_regulons <- GetTFRegulons(seurat_obj, wgcna_name)
    tf_net <- GetTFNetwork(seurat_obj, wgcna_name)

    if(is.null(tf_regulons) || nrow(tf_regulons) == 0){
        stop("No regulon data found. Check the WGCNA results in the Seurat object.")
    }
    if(is.null(tf_net) || nrow(tf_net) == 0){
        stop("No TF network data found. Ensure ConstructTFNetwork has been run.")
    }

    # get module info
    modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 'grey')
    module_tfs <- subset(modules, gene_name %in% unique(tf_regulons$tf))
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']

    # subset tf_regulons by module genes
    tf_regulons <- subset(tf_regulons, tf %in% modules$gene_name & gene %in% modules$gene_name)

    # use all genes or only TFs?
    if(TFs_only){
        tf_regulons <- subset(tf_regulons, gene %in% unique(tf_regulons$tf))
    }  

    # add module info to the TF and target gene
    ix <- match(tf_regulons$tf, modules$gene_name)
    tf_regulons$source_module <- modules$module[ix]
    ix <- match(tf_regulons$gene, modules$gene_name)
    tf_regulons$target_module <- modules$module[ix]
    tf_regulons <- tf_regulons %>% 
        dplyr::rename(source=tf, target=gene) %>%
        mutate(Gain = Gain * sign(Cor))
    cur_network <- tf_regulons

    # make empty matrices to store the number of links between mods
    pos_mat <- matrix(0, length(mods), length(mods))
    rownames(pos_mat) <- mods; colnames(pos_mat) <- mods
    neg_mat <- matrix(0, length(mods), length(mods))
    rownames(neg_mat) <- mods; colnames(neg_mat) <- mods

    # Identify the number of pos & neg regulatory links between mods
    reg_df <- data.frame()
    combos <- expand.grid(mods, mods)
    for(i in 1:nrow(combos)){
        m1 <- as.character(combos[i,'Var1'])
        m2 <- as.character(combos[i, 'Var2'])

        # how many positive connections from m1 to m2:
        cur_pos <- subset(cur_network, target_module == m1 & source_module == m2 & Gain >= 0)
        cur_neg <- subset(cur_network, target_module == m1 & source_module == m2 & Gain < 0)

        cur_df <- data.frame(
            source = m1, target = m2,
            n_pos = nrow(cur_pos),
            sum_pos = sum(cur_pos$Gain),
            n_neg = nrow(cur_neg),
            sum_neg = sum(cur_neg$Gain)
        )
        reg_df <- rbind(reg_df, cur_df)
    }

    # calculate averages
    reg_df$mean_pos <- reg_df$sum_pos / reg_df$n_pos
    reg_df$mean_pos <- ifelse(is.na(reg_df$mean_pos), 0, reg_df$mean_pos)
    reg_df$mean_neg <- reg_df$sum_neg / reg_df$n_neg
    reg_df$mean_neg <- ifelse(is.na(reg_df$mean_neg), 0, reg_df$mean_neg)

    # normalize counts by the total number of TFs in each module 
    tmp <- table(module_tfs$module)
    ix <- match(reg_df$source, names(tmp))
    reg_df$score_pos <- reg_df$n_pos / as.numeric(tmp)[ix]
    reg_df$score_pos <- ifelse(is.na(reg_df$score_pos), 0, reg_df$score_pos)
    reg_df$score_neg <- reg_df$n_neg / as.numeric(tmp)[ix]
    reg_df$score_neg <- ifelse(is.na(reg_df$score_neg), 0, reg_df$score_neg)

    # set factor levels for modules
    reg_df$source <- factor(reg_df$source, levels=mods)
    reg_df$target <- factor(reg_df$target, levels=mods)

    # return 
    reg_df

}


#' ModuleRegulatoryNetworkPlot
#'
#' This function visualizes the regulatory network between gene modules from a Seurat object. 
#' The plot displays transcription factor (TF) regulatory interactions between modules, with edges 
#' representing positive or negative regulatory links.
#'
#' @param seurat_obj A Seurat object containing single-cell data and WGCNA results.
#' @param feature A character string specifying the type of regulatory score to plot. Options are 
#' 'positive' (positive regulatory score), 'negative' (negative regulatory score), or 'delta' 
#' (difference between positive and negative scores). Default is 'delta'.
#' @param TFs_only Logical; if TRUE (default), only transcription factor (TF) genes are included 
#' in the plot. If FALSE, all genes are considered in the regulatory network.
#' @param layout A character string specifying the layout for the plot. Default is 'umap', which 
#' arranges the modules according to UMAP coordinates. Other layout options from \code{ggraph} 
#' can also be used.
#' @param umap_background Logical; if TRUE, the UMAP of module eigengenes is plotted in the 
#' background (if \code{layout = 'umap'}). Default is FALSE.
#' @param max_val Numeric; sets the maximum absolute value for the regulatory score. Any values 
#' exceeding this threshold are capped. Default is 1.
#' @param cutoff Numeric; edges with absolute regulatory scores below this value are excluded 
#' from the plot. Default is 0.
#' @param focus_source Character vector; optionally restricts the plot to only show regulatory 
#' links from specified source modules. Default is NULL (all modules included).
#' @param focus_target Character vector; optionally restricts the plot to only show regulatory 
#' links targeting specified modules. Default is NULL (all modules included).
#' @param loops Logical; if TRUE (default), loops (self-regulatory connections) are shown in the plot.
#' @param label_modules Logical; if TRUE (default), module names are displayed as labels in the plot.
#' @param high_color Character string; the color representing high regulatory scores in the 
#' color gradient. Default is 'orange2'.
#' @param mid_color Character string; the color representing intermediate regulatory scores in the 
#' color gradient. Default is 'white'.
#' @param low_color Character string; the color representing low regulatory scores in the 
#' color gradient. Default is 'dodgerblue'.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Additional arguments passed to the \code{ggraph::create_layout} function.
#'
#' @details
#' The function visualizes the regulatory network between gene modules by plotting transcription factor 
#' regulatory interactions. Positive and negative regulatory scores are displayed as edges connecting 
#' nodes (modules), with node sizes proportional to the number of genes in each module. Users can customize 
#' the network layout (e.g., UMAP or other graph layouts), filter edges based on score thresholds, and highlight 
#' specific modules by using the focus parameters. Edge colors represent the strength and direction (positive 
#' or negative) of regulatory interactions, and loops can optionally be plotted to indicate self-regulation.
#'
#' @return A ggplot2 object visualizing the regulatory network, which can be further customized or displayed 
#' using \code{plot()}.
#'
#' @import Seurat 
#' @import igraph 
#' @import ggraph 
#' @import tidygraph
#' @export
ModuleRegulatoryNetworkPlot <- function(
    seurat_obj,
    feature = 'delta',
    TFs_only=TRUE,
    layout = 'umap',
    umap_background=FALSE,
    max_val=1,
    cutoff=0,
    focus_source=NULL,
    focus_target=NULL,
    loops = TRUE,
    label_modules=TRUE,
    high_color = 'orange2',
    mid_color = 'white',
    low_color = 'dodgerblue',
    wgcna_name=NULL,
    ...
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # check that feature is valid 
    if(! feature %in% c('positive', 'negative', 'delta')){
        stop('Invalid selection for feature. Valid choices are positive, negative, or delta.')
    }

    # Check focus_source input
    if(!is.null(focus_source) && !all(focus_source %in% mods)){
        stop("Invalid modules in focus_source. Ensure they match the modules in the WGCNA experiment.")
    }

    # Check focus_target input
    if(!is.null(focus_target) && !all(focus_target %in% mods)){
        stop("Invalid modules in focus_target. Ensure they match the modules in the WGCNA experiment.")
    }

    # Check if UMAP layout is available
    if(layout == 'umap'){
        umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
        if(is.null(umap_df)){
            stop("UMAP not found in seurat_obj. Please run RunModuleUMAP first.")
        }
    }

    # get the module regulatory network:
    plot_df <- ModuleRegulatoryNetwork(seurat_obj, TFs_only=TFs_only, wgcna_name=wgcna_name)

    # get the modules table
    modules <- GetModules(seurat_obj, wgcna_name) %>% 
        subset(module != 'grey') %>% dplyr::mutate(module = droplevels(module))

    # get module color scheme
    mods <- levels(modules$module); mods <- mods[mods != 'grey']
    mod_colors <- modules %>% dplyr::select(c(module, color)) %>%
        dplyr::distinct() %>% dplyr::arrange(module) 
    cp <- mod_colors$color; names(cp) <- as.character(mod_colors$module)
    n_genes <- table(modules$module)[mods]
    size_limits <- range(as.numeric(n_genes))

    # what feature are we using?
    if(feature == 'positive'){
        plot_df$score <- plot_df$score_pos 
        label <- 'Positive \nRegulatory \nScore'
    } else if(feature == 'negative'){
        plot_df$score <- plot_df$score_neg
        label <- 'Negative \nRegulatory \nScore'
    } else{
        plot_df$score <- plot_df$score_pos - plot_df$score_neg
        label <- 'Pos - Neg \nRegulatory \nScore'
    }

    # remove links below the cutoff
    cur_net <- plot_df %>% 
        subset(abs(score) >= cutoff) 

    if(!is.null(focus_source)){
        if(!all(focus_source %in% mods)){stop("Invalid selection for focus_source. Must select the names of modules in your hdWGCNA experiment.")}
        cur_net <- subset(cur_net, source %in% focus_source)
    }
    else if(!is.null(focus_target)){
        if(!all(focus_target %in% mods)){stop("Invalid selection for focus_target. Must select the names of modules in your hdWGCNA experiment.")}
        cur_net <- subset(cur_net, target %in% focus_target)
    }

    # arrange by the source module
    cur_net <- cur_net %>%
        dplyr::arrange(source)

    # restrict range of values for plotting
    cur_net$score <- ifelse(abs(cur_net$score) > max_val, max_val * sign(cur_net$score), cur_net$score)
    cur_net$score <- ifelse(is.na(cur_net$score), 0, cur_net$score)

    # remove any modules that aren't here:
    mods_keep <- unique(c(unique(as.character(cur_net$source)), unique(as.character(cur_net$target))))
    cur_net$source <- factor(as.character(cur_net$source), levels=mods_keep)
    cur_net$target <- factor(as.character(cur_net$target), levels=mods_keep)
    n_genes <- n_genes[mods_keep] # changed this 

    # make a tidygraph object
    graph <- tidygraph::as_tbl_graph(cur_net) %>% 
        tidygraph::activate(nodes) 

    # set up umap layout
    if(layout == 'umap'){
        centroid_df <- umap_df %>% 
            dplyr::group_by(module) %>%
            dplyr::summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)) %>% 
            subset(module %in% mods_keep)
        centroid_df$module <- factor(as.character(centroid_df$module), levels=mods_keep)
        centroid_df <- centroid_df %>% dplyr::arrange(module)

        # create layout (UMAP)
        umap_layout <- centroid_df %>% 
            dplyr::rename(c(x=UMAP1, y = UMAP2, name=module)) %>% 
            as.data.frame()
        rownames(umap_layout) <- 1:nrow(umap_layout)

        # create the layout
        lay <- ggraph::create_layout(graph, umap_layout)
    } else{
        V(graph)$name <- mods_keep
        lay <- ggraph::create_layout(graph, layout=layout, ...)
    }

    # add the number of genes per module as the size
    lay$n_genes <- as.numeric(n_genes)

    # 1: initialize ggraph 
    p <- ggraph(lay) 

    # plot the module UMAP underneath?
    if(layout == 'umap' & umap_background){
        p <- p + geom_point(
            inherit.aes=FALSE, data=umap_df, 
            aes(x=UMAP1, y=UMAP2, size=kME), 
            color=umap_df$color, alpha=0.3
        )
    }

    # 2: Network edges
    p <- p + geom_edge_fan(
        aes(color=score, alpha=abs(score)),
        arrow = arrow(length = unit(2, 'mm'), type='closed'), 
        end_cap = circle(3, 'mm')
    ) 

    # 2.1: loops
    if(loops){
        p <- p + geom_edge_loop(
            aes(color=score, alpha=abs(score)),
        )
    }

    # 3: Network nodes (modules)
    p <- p + geom_node_point(
        aes(fill=name, size=n_genes), shape=21, color='black'
    ) 

    # 4: Add module labels
    if(label_modules){
        p <- p + geom_node_label(
            aes(label=name), repel=TRUE, max.overlaps=Inf, 
            color='black', size=3
        ) 
    }

    # 5: color scheme, edit theme
    p <- p +  scale_edge_colour_gradient2(high=high_color, mid=mid_color, low=low_color)  + 
    scale_colour_manual(values=cp) + 
    scale_fill_manual(values=cp) + 
    scale_size(limits=c(0, max(size_limits))) +
    labs(edge_colour=label) +
     guides(
        edge_alpha="none", 
        size = "none",
        fill = "none"
    ) 

    p

}

#' ModuleRegulatoryHeatmap
#'
#' This function visualizes the regulatory network between gene modules as a heatmap, where 
#' each cell represents the regulatory score between a source and target module. The heatmap 
#' can display either positive, negative, or delta (positive minus negative) regulatory scores. 
#' Optionally, a dendrogram can be plotted to cluster modules based on their regulatory interactions.
#'
#' @param seurat_obj A Seurat object containing single-cell data and WGCNA results.
#' @param feature A character string specifying the type of regulatory score to plot. Options are 
#' 'positive' (positive regulatory score), 'negative' (negative regulatory score), or 'delta' 
#' (difference between positive and negative scores). Default is 'delta'.
#' @param TFs_only Logical; if TRUE (default), only transcription factor (TF) genes are included 
#' in the heatmap. If FALSE, all genes are considered in the regulatory network.
#' @param dendrogram Logical; if TRUE (default), a dendrogram is added to the heatmap to cluster 
#' modules based on their regulatory interactions. Default is TRUE.
#' @param coord_fixed Logical; if TRUE (default), aspect ratio in x and y axes are equal. Default is TRUE.
#' @param max_val Numeric; sets the maximum absolute value for the regulatory score. Any values 
#' exceeding this threshold are capped. Default is 1.
#' @param min_val_label Numeric; the minimum number of interactions required for a label to be 
#' displayed on a heatmap cell. Default is 3.
#' @param high_color Character string; the color representing high regulatory scores in the 
#' heatmap. Default is 'orange2'.
#' @param mid_color Character string; the color representing intermediate regulatory scores in the 
#' heatmap. Default is 'white'.
#' @param low_color Character string; the color representing low regulatory scores in the 
#' heatmap. Default is 'dodgerblue'.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' 
#' @details
#' The function visualizes the regulatory network between gene modules by plotting regulatory 
#' scores as a heatmap. Each cell in the heatmap represents the regulatory interaction between a 
#' source module (columns) and a target module (rows). The color of the cell represents the magnitude 
#' and direction (positive or negative) of the regulatory score. Users can choose to plot positive, 
#' negative, or delta scores and can adjust the color gradient and score thresholds. A dendrogram 
#' can be added to cluster modules based on their regulatory patterns, and labels can be shown for 
#' cells with a sufficient number of interactions.
#'
#' @return A ggplot2 object visualizing the module regulatory network as a heatmap, which can be 
#' further customized or displayed using \code{plot()}.
#'
#' @import Seurat
#' @import ggraph
#' @import tidygraph
#' @export
ModuleRegulatoryHeatmap <- function(
    seurat_obj,
    feature = 'delta',
    TFs_only=TRUE,
    dendrogram=TRUE,
    coord_fixed = TRUE,
    max_val=1,
    min_val_label=3,
    high_color = 'orange2',
    mid_color = 'white',
    low_color = 'dodgerblue',
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # check that feature is valid 
    if(! feature %in% c('positive', 'negative', 'delta')){
        stop('Invalid selection for feature. Valid choices are positive, negative, or delta.')
    }

    # get the module regulatory network:
    plot_df <- ModuleRegulatoryNetwork(seurat_obj, TFs_only=TFs_only, wgcna_name=wgcna_name)

    # what feature are we using?
    if(feature == 'positive'){
        plot_df$score <- plot_df$score_pos 
        plot_df$n <- plot_df$n_pos
        label <- 'Positive \nRegulatory \nScore'
    } else if(feature == 'negative'){
        plot_df$score <- plot_df$score_neg
        plot_df$n <- plot_df$n_neg  
        label <- 'Negative \nRegulatory \nScore'
    } else{
        plot_df$score <- plot_df$score_pos - plot_df$score_neg
        plot_df$n <- 0
        label <- 'Pos - Neg \nRegulatory \nScore'
    }

    # restrict range of values for plotting
    plot_df$score <- ifelse(abs(plot_df$score) > max_val, max_val * sign(plot_df$score), plot_df$score)
    plot_df$score <- ifelse(is.na(plot_df$score), 0, plot_df$score)

    # labels
    plot_df$label <- ifelse(plot_df$n >= min_val_label, plot_df$n, "")

    # plot the dendrogram?
    if(dendrogram){

        tmp <- plot_df %>% 
            dplyr::select(c(source, target, score)) %>%
            tidyr::pivot_wider(names_from=source, values_from=score) %>%
            as.data.frame()

        rownames(tmp) <- tmp$target 
        tmp <- as.matrix(tmp[,-1])

        # cluster 
        hc <- hclust(dist(tmp))

        # make the dendrogram plot
        p_dendro <- ggraph::ggraph(hc, "dendrogram", height=height) + 
            ggraph::geom_edge_elbow() +
            theme(plot.margin = margin(c(0,0,0,0)))

        # re-order modules by dendrogram
        dend_order <- hc$label[hc$order]
        plot_df$source <- factor(as.character(plot_df$source), levels=dend_order)
        plot_df$target <- factor(as.character(plot_df$target), levels=dend_order)

    }

    # plot heatmap
    p <- plot_df %>% 
        ggplot(aes(x=source, y=fct_rev(target), fill=score)) + 
        geom_tile() 
        
    # add labels
    if(feature != "delta"){
        p <- p + geom_text(aes(label=label)) +
            scale_fill_gradient(low=mid_color, high=high_color) 
    } else{
        p <- p + scale_fill_gradient2(low=low_color, mid=mid_color, high=high_color)
    }

    # theme etc 
    p <- p +
        xlab('Source Module') + ylab('Target Module') + 
        Seurat::RotatedAxis() +
        labs(fill=label) +
        theme(
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(linewidth=1, fill=NA, color='black')
        )

    if(coord_fixed){
        p <- p + coord_fixed() 
    }

    if(dendrogram){
        p <- (p_dendro / p) + patchwork::plot_layout(heights=c(0.2, 1))
    }
    p

}
