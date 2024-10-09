#' AssignTFRegulons
#'
#' Define the set of likely target genes (Regulons) for each transcrition factor
#'
#' @return seurat_obj with the TF Regulon information added
#'
#' @param seurat_obj A Seurat object
#' @param strategy method for defining regulons, "A", "B", or "C". See Details for more info.
#' @param reg_thresh threshold for regulatory score in strategies A, B, and C
#' @param n_tfs for strategy A, the number of top TFs to keep for each gene
#' @param n_genes for strategy B, the number of top target genes to keep for each TF
#' @param wgcna_name name of the WGCNA experiment
#' @details 
#' AssignTFRegulons uses the TF network information from ConstructTFNetwork to define 
#' sets of confident TF-gene pairs. A "regulon" is the set of target genes for a given TF,
#' and this function provides three different strategies to define TF regulons. Strategy "A"
#' selects the top TFs for each gene, strategy "B" selects the top genes for each TF, and 
#' strategy "C" retains all TF-gene pairs above a certain regulatory score (reg_thresh). 
#' 
#' @import Seurat
#' @import Matrix
#' @export
AssignTFRegulons <- function(
    seurat_obj,
    strategy = "A", # A, B, or C
    reg_thresh = 0.01,
    n_tfs = 10,
    n_genes = 50,
    wgcna_name=NULL
){

    # get data from active assay if wgcna_name is not given
    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # check strategy
    if(!strategy %in% c("A", "B", "C")){
        stop("Invalid strategy. Choose one of: 'A', 'B', or 'C'.")
    }

    # check reg_thresh
    if(!is.numeric(reg_thresh) || reg_thresh < 0){
        stop("The 'reg_thresh' parameter must be a non-negative numeric value.")
    }

    # get the tf_net from the seurat_obj 
    tf_net <- GetTFNetwork(seurat_obj, wgcna_name)

    # check that the TF network is not empty
    if(is.null(tf_net) || nrow(tf_net) == 0){
        stop("The TF network is empty or not found. Ensure 'ConstructTFNetwork' was run successfully.")
    }

    # assign regulons based on the selected strategy
    if(strategy == 'A'){

        # Take the top TFs for each gene
        tf_regulons <- tf_net %>% 
            subset(Gain >= reg_thresh) %>% 
            dplyr::group_by(gene) %>%
            dplyr::slice_max(order_by=Gain, n=n_tfs) %>% 
            dplyr::ungroup()

    } else if(strategy == 'B'){

        # Take the top target genes for each TF
        tf_regulons <- tf_net %>% 
            subset(Gain >= reg_thresh) %>% 
            dplyr::group_by(tf) %>%
            dplyr::slice_max(order_by=Gain, n=n_genes) %>% 
            dplyr::ungroup()

    } else if(strategy == 'C'){

        # Take all interactions above a certain score
        tf_regulons <- tf_net %>% 
            subset(Gain >= reg_thresh) 

    } 

    # arrange the resulting regulons:
    tf_regulons <- tf_regulons %>% 
        dplyr::group_by(tf) %>% 
        dplyr::arrange(desc(Gain*sign(Cor)), .by_group=TRUE)

    # add the regulons to the seurat object
    tf_regulons <- as.data.frame(tf_regulons)
    seurat_obj <- SetTFRegulons(seurat_obj, tf_regulons, wgcna_name)

}


#' RegulonScores
#'
#' Calculate expression scores for TF Regulons
#'
#' @return seurat_obj with the TF Regulon Scores added
#'
#' @param seurat_obj A Seurat object
#' @param target_type Which types of TF target genes to compute scores for? "positive", "negative", or "both"
#' @param cor_thresh threshold for TF-gene correlation for genes to be included in the regulon score
#' @param exclude_grey_genes option to exclude genes that are in the grey module from the regulon scores
#' @param wgcna_name name of the WGCNA experiment
#' @details 
#' RegulonScores calculates expression signatures for each TF regulon using the UCell algorithm
#' This function can calculate separate scores for TF regulons for target genes that are positively
#' or negatively correlated with the TF, representing putative acivated or repressed genes. These scores 
#' conveniently summarize the expression levels of the entire TF regulon, similar to module eigengenes
#' for the co-expression network analyssis. 
#' 
#' @import Seurat
#' @import Matrix 
#' @import UCell
#' @export
RegulonScores <- function(
    seurat_obj,
    target_type = 'positive',
    cor_thresh = 0.05,
    exclude_grey_genes = TRUE,
    wgcna_name = NULL,
    ... # options to pass to UCell
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # get modules:
    modules <- GetModules(seurat_obj, wgcna_name)

    # get TF target genes:
    tf_regulons <- GetTFRegulons(seurat_obj, wgcna_name)

    # subset by type of target gene?
    if(target_type == 'positive'){
        tf_regulons <- subset(tf_regulons, Cor > cor_thresh)
    } else if(target_type == 'negative'){
        tf_regulons <- subset(tf_regulons, Cor < cor_thresh)
    } else if(target_type == 'both'){
        tf_regulons <- subset(tf_regulons, abs(Cor) > cor_thresh)
    } else{
        stop("Invalid selection for target_type. Valid choices are positive, negative, both.")
    }

    if(exclude_grey_genes){

        # subset modules
        modules <- subset(modules, module != 'grey')

        # subset regulons
        tf_regulons <- subset(tf_regulons, gene %in% modules$gene_name & tf %in% modules$gene_name)

    }

    # set up the lists
    tfs_use <- unique(tf_regulons$tf)
    target_genes <- lapply(tfs_use, function(cur_tf){
        subset(tf_regulons, tf == cur_tf) %>% .$gene
    })
    names(target_genes) <- tfs_use

    # use UCell to comptue the TF regulons cores
    regulon_scores <- UCell::AddModuleScore_UCell(
        seurat_obj, features=target_genes,
        ...
    )@meta.data
    regulon_scores <- regulon_scores[,paste0(tfs_use, '_UCell')]

    # rename the columns to remove "_UCell"
    colnames(regulon_scores) <- gsub("_UCell", "", colnames(regulon_scores))

    # add the regulon scores to the seurat object
    seurat_obj <- SetRegulonScores(
        seurat_obj, 
        regulon_scores,
        target_type,
        wgcna_name
    )

    seurat_obj
}


#' GetTFTargetGenes
#'
#' Retrieve target genes for a set of transcription factors (TFs) from a Seurat object, 
#' based on either regulons or the full TF network. The function allows exploration 
#' of direct targets or extending the network by adding target genes of TFs that are 
#' downstream in the network.
#'
#' @param seurat_obj A Seurat object containing the TF network or regulon information. 
#' This must include WGCNA results in the @misc slot.
#' @param selected_tfs A list of transcription factors (TFs) to use as the starting point 
#' for network exploration. 
#' @param depth Integer value specifying the number of layers to extend the TF network 
#' from the selected TFs. For example, if depth=2, the target genes of selected_tfs are shown,
#' and additional target genes are shown for TFs that are targets of the original selected_tfs.
#' @param target_type The type of target genes to include in the output. Options are "positive" 
#' (genes with expression positively correlated with TFs), "negative" (genes with expression 
#' negatively correlated with TFs), or "both" (default, includes all targets).
#' @param use_regulons Logical flag indicating whether to use regulons (default = TRUE) or 
#' the entire TF network for target gene discovery.
#' @param wgcna_name Character string specifying the name of the WGCNA experiment. This should 
#' be stored in the @misc slot of the Seurat object.
#' 
#' @details 
#' This function retrieves direct and extended targets of a set of TFs based on the regulatory 
#' network stored in the Seurat object. The depth parameter controls how many layers of 
#' TF-target interactions to explore, while the target_type parameter allows filtering 
#' based on correlation direction. The use_regulons flag specifies whether to use regulon 
#' data, which may provide a more confident set of interactions, or the full TF network.
#' 
#' @return A data frame containing the TF-target interactions at each specified depth level, 
#' with additional information such as regulatory score, correlation, and depth.
#' 
#' @import dplyr
#' @import Seurat
#' @export
GetTFTargetGenes <- function(
    seurat_obj, 
    selected_tfs,
    depth=1,
    target_type = 'both', 
    use_regulons=TRUE,
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # check depth
    if(!is.numeric(depth) || depth < 1 || depth %% 1 != 0){
        stop("Invalid input: depth must be a positive integer.")
    }

    # check target types
    valid_target_types <- c("positive", "negative", "both")
    if(!(target_type %in% valid_target_types)){
        stop(paste0("Invalid target_type: must be one of ", paste(valid_target_types, collapse=', '), "."))
    }

    # get the regulons (recommended), or use the full TF network?
    if(use_regulons){
        tf_regulons <- GetTFRegulons(seurat_obj, wgcna_name)
    } else{
        tf_regulons <- GetTFNetwork(seurat_obj, wgcna_name)
    }

    # check if selectd_tfs is valid 
    if(! all(selected_tfs %in% tf_regulons$tf)){
        not_found <- selected_tfs[which(!selected_tfs %in% tf_regulons$tf)]
        stop(paste0("The following TFs were not found in the TF network: ",  paste0(not_found, collapse=' ,')))
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

    # initial condition
    prev_tfs <- selected_tfs
    cur_network <- data.frame()

    # loop for each depth level
    for(i in 1:depth){

        cur_regulons <- tf_regulons %>% 
            subset(tf %in% prev_tfs) %>% 
            mutate(depth = i)
        cur_network <- rbind(cur_network,  cur_regulons)

        # subset for tfs
        cur_tfs <- cur_regulons %>% 
            subset(gene %in% unique(tf_regulons$tf)) %>% .$gene

        prev_tfs <- unique(c(prev_tfs, cur_tfs))
    }
    
    return(cur_network)

}


#' FindDifferentialRegulons
#'
#' Function to compare TF regulon scores between two sets of cell barcodes.
#'
#' @return A dataframe contaning differential regulon results
#'
#' @param seurat_obj A Seurat object
#' @param barcodes1 character vector containing cell barcodes for the first group to test. Positive fold-change means up-regulated in this group.
#' @param barcodes2 character vector containing cell barcodes for the second group to test. Negative fold-change means up-regulated in this group.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param assay Assay to extract data for aggregation. Default = 'RNA'
#' @param slot Slot to extract data for aggregation. Default = 'counts'. Slot is used with Seurat v4 instead of layer.
#' @param layer Layer to extract data for aggregation. Default = 'counts'. Layer is used with Seurat v5 instead of slot.
#' @param ... Additional parameters for the Seurat FindMarkers function
#' 
#' @details
#' FindDifferentialRegulons compares two groups based on their TF regulon scores. 
#' Three comparisons are made for different sets of features: positive regulon scores, negative regulon scores, and gene expression.
#' The same settings for the test will be used for all three tests.
#' 
#' @export
FindDifferentialRegulons <- function(
    seurat_obj,
    barcodes1,
    barcodes2,
    assay = 'RNA',
    slot = 'data',
    layer = 'data',
    wgcna_name=NULL,
    test.use='wilcox',
    only.pos=FALSE,
    logfc.threshold = 0,
    min.pct=0,
    verbose=FALSE,
    pseudocount.use=0,
    ...
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # check regulon scores
    tmp1 <- GetRegulonScores(seurat_obj, 'positive', wgcna_name)
    tmp2 <- GetRegulonScores(seurat_obj, 'negative', wgcna_name)
    if(is.null(tmp1) | is.null(tmp2)){
        stop('Regulon scores not found, please run RegulonScores first for target_type="positive" and target_type="negative".' )
    }

    # ensure that selected barcodes are in the seurat obj
    if(!(all(barcodes1 %in% colnames(seurat_obj)))){
    stop('Some barcodes in barcodes1 not found in colnames(seurat_obj).')
    }
    if(!(all(barcodes2 %in% colnames(seurat_obj)))){
    stop('Some barcodes in barcodes2 not found in colnames(seurat_obj).')
    }

    # check for overlap in the two groups of bacrodes:
    if(length(intersect(barcodes1, barcodes2)) > 0){
    stop('Some barcodes overlap in barcodes1 and barcodes2')
    }

    dregs_list <- list()
    for(target_type in c('positive', 'negative')){

        # get regulon scores
        regulon_scores <- GetRegulonScores(
            seurat_obj, 
            target_type = target_type, 
            wgcna_name = wgcna_name
        )

        # create new assay for regulons:
        reg_assay <- Seurat::CreateAssayObject(t(regulon_scores))

        # run differential test on the reg assay
        cur_dregs <- FindMarkers(
            reg_assay,
            cells.1 = barcodes1,
            cells.2 = barcodes2,
            slot='counts', # should I make this layer at some point?
            test.use=test.use,
            only.pos=only.pos,
            logfc.threshold=logfc.threshold,
            min.pct=min.pct,
            verbose=verbose,
            pseudocount.use=pseudocount.use,
            ...
        )

        # rename:
        cur_dregs <- cur_dregs %>%
            dplyr::select(-c(pct.1, pct.2)) 
        colnames(cur_dregs) <- paste0(colnames(cur_dregs), '_', target_type) 

        cur_dregs$tf <- rownames(cur_dregs)

        dregs_list[[target_type]] <- cur_dregs

    }

    # merge:
    dregs <- dplyr::left_join(
        dregs_list[[1]],
        dregs_list[[2]],
        by='tf'
    )

    tf_genes <- dregs$tf

    # get expression matrix:
    if(CheckSeurat5()){
        X <- SeuratObject::LayerData(seurat_obj, assay=assay, layer=layer)
    } else{
        X <- Seurat::GetAssayData(seurat_obj, assay=assay, slot=slot)
    }
    X <- X[tf_genes,]
    exp_assay <- Seurat::CreateAssayObject(X)

    degs <- FindMarkers(
        exp_assay,
        cells.1 = group1,
        cells.2 = group2,
        slot='data',
        test.use=test.use,
        only.pos=only.pos,
        logfc.threshold=logfc.threshold,
        min.pct=min.pct,
        verbose=verbose,
        pseudocount.use=pseudocount.use,
        ...
    )
    degs$gene <- rownames(degs)

    # join them together!!!
    tmp <- dplyr::inner_join(
        dregs %>% dplyr::rename(gene = tf), 
        degs %>% dplyr::select(c(p_val_adj, avg_log2FC, gene)), 
        by = 'gene'
    ) 
    dregs$p_val_deg <- tmp$p_val
    dregs$avg_log2FC_deg <- tmp$avg_log2FC
    dregs$p_val_adj_deg <- tmp$p_val_adj 

    # move tf to first col
    dregs <- dregs %>% dplyr::relocate(tf)

    # add additional info
    hub_df <- GetHubGenes(seurat_obj, n_hubs=Inf, wgcna_name=wgcna_name)
    rownames(hub_df) <- hub_df$gene_name
    dregs$module <- hub_df[dregs$tf, 'module']
    dregs$kME <- hub_df[dregs$tf, 'kME']

    # return
    dregs

}

#' PlotDifferentialRegulons
#'
#' Function to visualize differential TF regulon activity between two groups of cells based on regulon scores.
#' The plot shows the average log fold-change (logFC) for positive and negative regulon scores in a scatter plot, 
#' with the points colored by their associated module and optionally labeled with the most differentially active regulons.
#'
#' @param seurat_obj A Seurat object that contains the WGCNA experiment and TF regulon scores.
#' @param dregs Dataframe output from the \code{FindDifferentialRegulons} function containing differential regulon results.
#' @param n_label Integer specifying the number of top up- and down-regulated regulons to label on the plot. If 'all', all significant regulons will be labeled. Default = 10.
#' @param logfc_thresh Numeric threshold for labeling and annotating regulons based on logFC values. Default = 0.1.
#' @param lm Logical indicating whether to plot a linear regression line for the logFC values of positive vs. negative regulons. Default = TRUE.
#' @param wgcna_name The name of the WGCNA experiment in the \code{seurat_obj@misc} slot. Default is the active WGCNA experiment.
#'
#' @details
#' This function generates a scatter plot where each point represents a TF regulon, with the x-axis showing the average log fold-change for positive regulon scores and the y-axis showing the average log fold-change for negative regulon scores. Points are colored by their module assignment, and the size of the points reflects the module eigengene correlation (kME) value of the TF.
#' 
#' Differentially expressed regulons can be highlighted, and the most significantly up- and down-regulated regulons can be labeled. Additional options allow adding a linear regression line and controlling label density based on logFC thresholds.
#'
#' @return A ggplot object representing the differential regulon scatter plot.
#'
#' @export
PlotDifferentialRegulons <- function(
    seurat_obj,
    dregs,
    n_label=10,
    logfc_thresh=0.1,
    lm = TRUE,
    wgcna_name = NULL
){

    # checks
    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    if(!is.data.frame(dregs)){
        stop("The input 'dregs' must be a dataframe, output from the FindDifferentialRegulons function.")
    }

    if(!is.numeric(n_label) || n_label < 1){
        if(n_label != "all"){
            stop("The 'n_label' parameter must be a positive integer or 'all'.")
        }
    }

    if(!is.numeric(logfc_thresh) || logfc_thresh < 0){
        stop("The 'logfc_thresh' parameter must be a non-negative numeric value.")
    }

    # get the module color scheme
    modules <- GetModules(seurat_obj, wgcna_name) %>% 
        subset(module != 'grey') %>% 
        dplyr::mutate(module=droplevels(module))
    module_colors <- modules %>% 
        dplyr::select(c(module, color)) %>% 
        dplyr::distinct()
    mods <- levels(modules$module)
    mods <- mods[mods %in% dregs$module]
    mod_colors <- module_colors$color; names(mod_colors) <- as.character(module_colors$module)

    # x and y axis limits
    plot_range <- max(c(abs(range(dregs$avg_log2FC_positive)), abs(range(dregs$avg_log2FC_negative))))

    # annotations for the corners 
    up_right <- dregs %>% subset(avg_log2FC_positive >= logfc_thresh & avg_log2FC_negative >= logfc_thresh & (p_val_adj_negative <= 0.05 | p_val_adj_positive <= 0.05)) %>% nrow
    down_right <- dregs %>% subset(avg_log2FC_positive >= logfc_thresh & avg_log2FC_negative <= -logfc_thresh & (p_val_adj_negative <= 0.05 | p_val_adj_positive <= 0.05)) %>% nrow
    up_left <- dregs %>% subset(avg_log2FC_positive <= -logfc_thresh & avg_log2FC_negative >= logfc_thresh & (p_val_adj_negative <= 0.05 | p_val_adj_positive <= 0.05)) %>% nrow
    down_left <- dregs %>% subset(avg_log2FC_positive <= -logfc_thresh & avg_log2FC_negative <= -logfc_thresh & (p_val_adj_negative <= 0.05 | p_val_adj_positive <= 0.05)) %>% nrow

    annotations <- data.frame(
            xpos = c(-Inf,-Inf,Inf,Inf),
            ypos =  c(-Inf, Inf,-Inf,Inf),
            annotateText = c(as.character(down_left),as.character(up_left), as.character(down_right),as.character(up_right)),
            hjustvar = c(-1,-1,2,2),
            vjustvar = c(-1,2,-1,2)) #<- adjust

    up_genes <- dregs %>% subset(avg_log2FC_positive <= -logfc_thresh & avg_log2FC_negative >= logfc_thresh & (p_val_adj_negative <= 0.05 | p_val_adj_positive <= 0.05)) %>% .$tf
    down_genes <- dregs %>% subset(avg_log2FC_positive >= logfc_thresh & avg_log2FC_negative <= -logfc_thresh & (p_val_adj_negative <= 0.05 | p_val_adj_positive <= 0.05)) %>% .$tf
    signif_genes <- c(up_genes, down_genes)
    if(n_label == 'all'){
        label_genes <- signif_genes
    } else{
        up_genes <- dregs %>% subset(avg_log2FC_positive <= -logfc_thresh & avg_log2FC_negative >= logfc_thresh & (p_val_adj_negative <= 0.05 | p_val_adj_positive <= 0.05)) %>% dplyr::slice_max(n=n_label, order_by=avg_log2FC_negative - avg_log2FC_positive) %>% .$tf
        down_genes <- dregs %>% subset(avg_log2FC_positive >= logfc_thresh & avg_log2FC_negative <= -logfc_thresh & (p_val_adj_negative <= 0.05 | p_val_adj_positive <= 0.05)) %>% dplyr::slice_max(n=n_label, order_by=avg_log2FC_positive - avg_log2FC_negative) %>% .$tf
        label_genes <- c(up_genes, down_genes)
    }
    
    dregs$label <- ifelse(dregs$tf %in% label_genes, dregs$tf, '')
    dregs$signif <- ifelse(dregs$tf %in% signif_genes, dregs$tf, '')

    # denote DEGs 
    degs <- dregs %>% subset(abs(avg_log2FC_deg) >= logfc_thresh & p_val_adj_deg < 0.05) %>% .$tf
    dregs$deg <- dregs$tf %in% degs

    # make the plot
    p <- dregs %>%
        ggplot(aes(x=avg_log2FC_positive, y=avg_log2FC_negative, color=module, fill=module, size=kME)) + 
        geom_hline(yintercept=0, linetype='dashed', color='black') +
        geom_vline(xintercept=0, linetype='dashed', color='black') 
    
    # add the unlabeled points
    p <- p + 
        geom_point(
            data = subset(dregs, signif == "" & !deg), 
            alpha=0.5) +
        geom_point(
            data = subset(dregs, signif == "" & deg), 
            alpha=0.5, shape=18) 

    # add the labeled points
    p <- p + geom_point(
            data = subset(dregs, signif != "" & !deg), 
            color='black', shape=21
        ) + geom_point(
        data = subset(dregs, signif != "" & deg), 
        color='black', shape=23
    ) 
    
    # linear regression line
    if(lm){
        p <- p + 
            geom_smooth(
                inherit.aes=FALSE, 
                data=dregs, 
                mapping = aes(x = avg_log2FC_positive, y = avg_log2FC_negative), 
                method='lm', color='black'
            ) 
    }

    # add labels & modify appearance
    p <- p +  
        geom_text_repel(aes(label=label), max.overlaps=Inf, color='black', size=3) +
        geom_text(inherit.aes=FALSE, data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
        scale_color_manual(values=mod_colors) +
        scale_fill_manual(values=mod_colors) +
        xlim(c(-plot_range, plot_range)) +
        ylim(c(-plot_range, plot_range)) +
        xlab(bquote("Avg. log"[2]~"(FC), Positive Regulons")) +
        ylab(bquote("Avg. log"[2]~"(FC), Negative Regulons")) +
        theme(
            panel.border = element_rect(color='black', fill=NA, size=1),
            panel.grid.major = element_blank(),
            axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position='bottom'
        ) + Seurat::NoLegend() 

    p

}


