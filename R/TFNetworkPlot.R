#' TFNetworkPlot
#'
#' Plots a network of TFs and predicted target genes.
#'
#' @return ggplot object containing the TFNetworkPlot
#'
#' @param seurat_obj A Seurat object
#' @param selected_tfs A list of TFs 
#' @param depth Number of layers to extend the TF network from the selected_tfs. For example, if depth=2 (default), the target genes of selected_tfs are shown, and additional target genes are shown for other TFs that are target genes of the original selected_tfs.
#' @param edge_weight Attribute to use to color the network edges. "Cor" to show the pearson correlation, or "Gain" to show the importance score from the XGBoost model.
#' @param cutoff Cutoff for the edge weights, links below this value will not be included.
#' @param color_cutoff Maximum value for the colorscale for the edge weights.
#' @param target_type Type of target genes to show in the network. "Positive" shows genes with expression positively correlated with selected_tfs, "Negative" shows genes that are negatively correlated, and "Both" shows both (default).
#' @param use_regulons Logical indicating whether to use the TF-gene links defined in the TF Regulons (TRUE, default), or all putative TF-gene links (FALSE).
#' @param label_genes List of target genes to label in the network plot.
#' @param label_TFs The depth level in the network to plot the names of TFs. 
#' @param no_labels Logical indicating whether or not to remove all labels.
#' @param TFs_only Logical indicating whether or not to use only TFs in the plot or to use TFs and other genes.
#' @param high_color Color corresponding to positive edge weights
#' @param mid_color Color corresponding to edge weights close to 0
#' @param low_color Color corresponding to negative edge weights
#' @param node_colors List of color names for the TFs and genes at each depth. For example, if depth=2, you should supply a list of 3 colors corresponding to the selected_tfs (depth 0), the primary target genes (depth 1), and the secondary target genes (depth 2).
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' TFNetworkPlot plots a network of TFs and their predicted target genes, based on the results 
#' from ConstructTFNetwork.
#' 
#' @import Seurat
#' @import igraph
#' @import ggraph
#' @import tidygraph
#' @export
TFNetworkPlot <- function(
    seurat_obj,
    selected_tfs,
    depth = 2, 
    edge_weight = 'Cor', # you can use Cor or Gain
    cutoff = 0.01,
    color_cutoff = 0.75,
    target_type = 'both', # 'negative', 'all'
    use_regulons=TRUE,
    label_genes = NULL,
    label_TFs = 1, # label TFs at this depth and below
    no_labels = FALSE,
    TFs_only = FALSE,
    high_color = 'orange2',
    mid_color = 'white',
    low_color = 'dodgerblue',
    node_colors =  c('black', 'darkorchid4', 'mediumpurple2'),
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # do we have enough colors?
    if(length(node_colors) < depth){
        stop('Insufficient number of colors supplied in node_colors, should be equal to depth.')
    }

    # get the regulons (recommended), or use the full TF network?
    tf_net <- GetTFNetwork(seurat_obj, wgcna_name)
    if(use_regulons){
        tf_regulons <- GetTFRegulons(seurat_obj, wgcna_name)
    } else{
        tf_regulons <- tf_net
    }

    # check if edge_weight is valid
    if(! edge_weight %in% c('Cor', 'Gain')){
        stop('Invalid selection for edge_weight. Valid options are Cor or Gain.')
    }

    # check if selectd_tfs is valid 
    if(! all(selected_tfs %in% tf_regulons$tf)){
        not_found <- selected_tfs[which(!selected_tfs %in% tf_regulons$tf)]
        stop(paste0("The following TFs were not found in the TF network: ",  paste0(not_found, collapse=' ,')))
    }

    # check if label_genes is valid
    if(is.null(label_genes)){
        label_genes <- c()
    } else{
        if(!any(label_genes %in% tf_regulons$gene)){
            not_found <- label_genes[which(!label_genes %in% tf_regulons$gene)]
            warning("Some genes in label_genes not found in the TF network.")
            label_genes <- label_genes[label_genes %in% tf_regulons$gene]
        }
    }

    # subset by type of target gene?
    if(target_type == 'positive'){
        tf_regulons <- subset(tf_regulons, Cor > 0)
    } else if(target_type == 'negative'){
        tf_regulons <- subset(tf_regulons, Cor < 0)
    } else if(target_type == 'both'){
        tf_regulons <- tf_regulons
    } else{
        stop("Invalid selection for target_type. Valid choices are positive, negative, both.")
    }

    # get the target genes
    cur_network <- GetTFTargetGenes(
        seurat_obj,
        selected_tfs=selected_tfs, 
        depth=depth, 
        target_type=target_type,
        use_regulons=use_regulons,
        wgcna_name=wgcna_name
    )

    # get the max depth of each gene
    gene_depths <- cur_network %>% 
        group_by(gene) %>% 
        slice_min(n=1, order_by=depth) %>% 
        select(c(gene, depth)) %>% distinct()

    # set the column that will be used as the edge weight
    cur_network$edge_weight <- cur_network[,edge_weight]
    if(edge_weight == 'Gain'){
        cur_network$edge_weight <- cur_network$edge_weight * sign(cur_network$Cor)
    }

    # cutoff for edge_weight
    cur_network$edge_weight <- ifelse(abs(cur_network$edge_weight) > color_cutoff, sign(cur_network$edge_weight) * color_cutoff, cur_network$edge_weight)
    cur_network <- subset(cur_network, abs(edge_weight) >= cutoff)

    # make an igraph network 
    cur_network <- cur_network %>%
        dplyr::rename(c(source=tf, target=gene)) %>%
        mutate(Score = sign(Cor) * Gain)

    # get a list of genes in the network    
    cur_genes <- unique(cur_network$target)

    # do we only want to plot TFs?
    if(TFs_only){cur_network <- subset(cur_network, target %in% unique(tf_net$tf))}

    # make a tidygraph
    graph <- tidygraph::as_tbl_graph(cur_network) %>% 
        tidygraph::activate(nodes) %>% 
        mutate(degree  = centrality_degree())  

    # compute the degree for each TF:
    tf_degrees <- table(tf_regulons$tf)
    tmp <- tf_degrees[names(V(graph))]; tmp <- tmp[!is.na(tmp)]
    V(graph)[names(tmp)]$degree <- as.numeric(tmp)

    # specify the selected TFs vs TFs vs genes
    V(graph)$gene_type <- ifelse(names(V(graph)) %in% unique(tf_regulons$tf), 'TF', 'Gene')
    V(graph)$gene_type <- ifelse(names(V(graph)) %in% selected_tfs, 'selected', V(graph)$gene_type)

    # set up the graph layout
    if(igraph::is_connected(graph)){
        n_pivots <- 250
        if(length(V(graph)) < n_pivots){
            n_pivots <- length(V(graph)) / 2
        }
        lay <- ggraph::create_layout(graph, layout='sparse_stress', pivots=n_pivots)
    } else{
        lay <- ggraph::create_layout(graph, layout='igraph', algorithm='nicely')
    }

    # correct for the edge case where there are Nas
    lay[is.na(lay$x),'x'] <- 0
    lay[is.na(lay$y),'y'] <- 0
    
    # add extra info
    lay$size <- ifelse(lay$name %in% unique(tf_regulons$tf), 5, 2)

    # add the depth:
    tmp <- dplyr::left_join(lay, gene_depths, by = c('name' = 'gene'))
    lay$depth <- tmp$depth
    lay$depth <- ifelse(lay$name %in% selected_tfs, 0, lay$depth)
    lay$depth <- factor(lay$depth, levels=0:depth)

    # shape layout:
    cur_shapes <- c(18, 17, 16); names(cur_shapes) <- c('selected', 'TF', 'Gene')

    # set up plotting labels
    label_tfs <- subset(cur_network, depth <= label_TFs & target %in% tf_regulons$tf) %>% .$target %>% unique
    lay$lab <- ifelse(lay$name %in% c(selected_tfs, label_tfs, label_genes), lay$name, NA)

    # node colors
    cp <- node_colors[1:(depth+1)]
    names(cp) <- levels(lay$depth)

    # color by depth
    p <- ggraph(lay) + 
        geom_edge_link(
            aes(color=edge_weight, alpha=abs(edge_weight)),
            arrow = arrow(length = unit(1, 'mm'), type='closed')
        ) + 
        geom_node_point(data=subset(lay, (! name %in% tf_regulons$tf)), aes(color=depth, shape=gene_type, size=degree)) +
        geom_node_point(data=subset(lay, name %in% tf_regulons$tf & !(name %in% selected_tfs)), aes(fill=depth, size=degree), color='black', shape=25) +
        geom_node_point(data=subset(lay, name %in% selected_tfs ), aes(color=depth, shape=gene_type, size=degree)) +
        geom_node_point(data=subset(lay, name %in% selected_tfs ), aes(fill=depth, size=degree), color = 'black', shape=23) +
        scale_edge_colour_gradient2(high=high_color, mid=mid_color, low=low_color)  + 
        scale_colour_manual(values=cp) + 
        scale_fill_manual(values=cp) + 
        scale_shape_manual(values=cur_shapes) 

    # add the labels
    if(! no_labels){
        p <- p +         
            geom_node_label(aes(label=lab), repel=TRUE, max.overlaps=Inf, fontface='italic')
    }

    # clean up the legends
    p <- p + guides(
        edge_alpha="none", 
        size = "none",
        shape = "none",
        fill = "none"
    ) 
    p <- p + labs(edge_colour='strength')

    p

}

