#' ConstructTFNetwork
#'
#' Construct a network of transcription factors and target genes based on gene co-expression
#'
#' @return seurat_obj with the TF network added
#'
#' @param seurat_obj A Seurat object
#' @param model_params a list of model parameters to pass to xgboost
#' @param nfold number of folds for cross validation
#' @param callbacks callback functions for the xgboost model
#' @param wgcna_name name of the WGCNA experiment
#'
#' @details
#' ConstructTFNetwork uses motif-gene information to build a directed network of transcription
#' factors (TFs) and target genes. XGBoost regression is leveraged to model the expression of
#' each gene based on its candidate TF regulators. This analysis gives us information about
#' how important each TF is for predicting each gene, allowing us to prioritize the most likely
#' regulators of each gene. This process is done on the gene expression matrix stored with SetDatExpr,
#' which is typically the hdWGCNA metacell gene expression matrix.
#'
#' @export
ConstructTFNetwork <- function(
  seurat_obj,
  model_params,
  nfold = 5,
  callbacks = NULL,
  wgcna_name = NULL
) {
  check_suggests('xgboost', feature = 'ConstructTFNetwork')

  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get the motif information from the Seurat object:
  motif_matrix <- GetMotifMatrix(seurat_obj)
  motif_df <- GetMotifs(seurat_obj)

  if (is.null(motif_df)) {
    stop('Motif info not found in the seruat_obj. Please run MotifScan first.')
  }

  # check that gene names column is in the motif names
  if (!'gene_name' %in% colnames(motif_df)) {
    stop('gene_name column missing in motif table (GetMotifs(seurat_obj)). Please add a column indicating the gene_name in the seurat_obj for each motif.')
  }

  check_xgboost1 <- startsWith(as.character(packageVersion('xgboost')), '1')

  # define XGBoost callback function based on the installed version, because after v1 they changed some parameters
  if (is.null(callbacks)) {
    # looked up by name: the function only exists in the installed xgboost version
    cb_name <- ifelse(check_xgboost1, 'cb.cv.predict', 'xgb.cb.cv.predict')
    callbacks <- list(getExportedValue('xgboost', cb_name)(save_models = TRUE))
  }

  # subset the motif_df by genes that are in the seurat obj:
  motif_df <- subset(motif_df, gene_name %in% rownames(seurat_obj))

  # get the expression matrix:
  datExpr <- as.matrix(GetDatExpr(seurat_obj, wgcna_name = wgcna_name))
  genes_use <- colnames(datExpr)
  genes_use <- genes_use[genes_use %in% rownames(motif_matrix)]

  # set up output dataframes
  importance_df <- data.frame()
  eval_df <- data.frame()

  # set up the progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(genes_use), style = 3, width = 50, char = '=')
  counter <- 1

  for (cur_gene in genes_use) {
    setTxtProgressBar(pb, counter)

    # get the list of TFs that regulate this gene:
    cur_tfs <- names(which(motif_matrix[cur_gene, ]))
    cur_tfs <- subset(motif_df, motif_ID %in% cur_tfs) %>%
      .$gene_name %>%
      unique()
    cur_tfs <- cur_tfs[cur_tfs %in% genes_use]

    # set up the expression matrices
    if (cur_gene %in% cur_tfs) {
      cur_tfs <- cur_tfs[cur_tfs != cur_gene]
      x_vars <- datExpr[, cur_tfs]
    }
    x_vars <- datExpr[, cur_tfs]
    y_var <- as.numeric(datExpr[, cur_gene])

    if (length(cur_tfs) < 2) {
      print(paste0('Not enough putative TFs, skipping ', cur_gene))
      next
    }

    # correlation:
    tf_cor <- as.numeric(cor(
      x = as.matrix(x_vars),
      y = y_var
    ))
    names(tf_cor) <- cur_tfs

    # run xgboost model
    if (all(y_var == 0)) {
      print(paste0('skipping ', cur_gene))
      next
    }

    xgb <- xgboost::xgb.cv(
      params = model_params,
      data = xgboost::xgb.DMatrix(x_vars, label = y_var),
      nrounds = 100,
      showsd = FALSE,
      nfold = nfold,
      callbacks = callbacks,
      verbose = FALSE
    )

    # get the CV evaluation info
    xgb_eval <- as.data.frame(xgb$evaluation_log)
    xgb_eval$variable <- cur_gene

    # average the importance score from each fold
    importance <- Reduce('+', lapply(1:nfold, function(i) {
      if (check_xgboost1) {
        cur_model <- xgb$models[[i]]
      } else {
        cur_model <- xgb$cv_predict$models[[i]]
      }
      cur_imp <- xgboost::xgb.importance(feature_names = colnames(x_vars), model = cur_model)
      ix <- match(colnames(x_vars), as.character(cur_imp$Feature))
      cur_imp <- as.matrix(cur_imp[ix, -1])
      cur_imp[is.na(cur_imp)] <- 0
      cur_imp
    })) / nfold
    importance <- as.data.frame(importance)

    # add tf and source info
    importance$tf <- colnames(x_vars)
    importance$gene <- cur_gene

    # add the tf correlation information
    importance$Cor <- as.numeric(tf_cor)

    # re-order columns, and re-order rows by gain:
    importance <- importance %>% dplyr::select(c(tf, gene, Gain, Cover, Frequency, Cor))
    importance <- arrange(importance, -Gain)

    # append
    importance_df <- rbind(importance_df, importance)
    eval_df <- rbind(eval_df, xgb_eval)

    # update progress bar
    counter <- counter + 1
  }

  # close the progress bar
  close(pb)

  # add the results to the Seurat object
  seurat_obj <- SetTFNetwork(seurat_obj, importance_df, wgcna_name = wgcna_name)
  seurat_obj <- SetTFEval(seurat_obj, eval_df, wgcna_name = wgcna_name)

  # return(list(importance=importance_df, eval=eval_df))
  seurat_obj
}


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
#' @export
TFNetworkPlot <- function(
  seurat_obj,
  selected_tfs,
  depth = 2,
  edge_weight = 'Cor', # you can use Cor or Gain
  cutoff = 0.01,
  color_cutoff = 0.75,
  target_type = 'both', # 'negative', 'all'
  use_regulons = TRUE,
  label_genes = NULL,
  label_TFs = 1, # label TFs at this depth and below
  no_labels = FALSE,
  TFs_only = FALSE,
  high_color = 'orange2',
  mid_color = 'white',
  low_color = 'dodgerblue',
  node_colors = c('black', 'darkorchid4', 'mediumpurple2'),
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # do we have enough colors?
  if (length(node_colors) < depth) {
    stop('Insufficient number of colors supplied in node_colors, should be equal to depth.')
  }

  # get the regulons (recommended), or use the full TF network?
  tf_net <- GetTFNetwork(seurat_obj, wgcna_name)
  if (use_regulons) {
    tf_regulons <- GetTFRegulons(seurat_obj, wgcna_name)
  } else {
    tf_regulons <- tf_net
  }

  # check if edge_weight is valid
  if (!edge_weight %in% c('Cor', 'Gain')) {
    stop('Invalid selection for edge_weight. Valid options are Cor or Gain.')
  }

  # check if selectd_tfs is valid
  if (!all(selected_tfs %in% tf_regulons$tf)) {
    not_found <- selected_tfs[which(!selected_tfs %in% tf_regulons$tf)]
    stop(paste0('The following TFs were not found in the TF network: ', paste0(not_found, collapse = ' ,')))
  }

  # check if label_genes is valid
  if (is.null(label_genes)) {
    label_genes <- c()
  } else {
    if (!any(label_genes %in% tf_regulons$gene)) {
      not_found <- label_genes[which(!label_genes %in% tf_regulons$gene)]
      warning('Some genes in label_genes not found in the TF network.')
      label_genes <- label_genes[label_genes %in% tf_regulons$gene]
    }
  }

  # subset by type of target gene?
  if (target_type == 'positive') {
    tf_regulons <- subset(tf_regulons, Cor > 0)
  } else if (target_type == 'negative') {
    tf_regulons <- subset(tf_regulons, Cor < 0)
  } else if (target_type == 'both') {
    tf_regulons <- tf_regulons
  } else {
    stop('Invalid selection for target_type. Valid choices are positive, negative, both.')
  }

  # get the target genes
  cur_network <- GetTFTargetGenes(
    seurat_obj,
    selected_tfs = selected_tfs,
    depth = depth,
    target_type = target_type,
    use_regulons = use_regulons,
    wgcna_name = wgcna_name
  )

  # get the max depth of each gene
  gene_depths <- cur_network %>%
    group_by(gene) %>%
    slice_min(n = 1, order_by = depth) %>%
    select(c(gene, depth)) %>%
    distinct()

  # set the column that will be used as the edge weight
  cur_network$edge_weight <- cur_network[, edge_weight]
  if (edge_weight == 'Gain') {
    cur_network$edge_weight <- cur_network$edge_weight * sign(cur_network$Cor)
  }

  # cutoff for edge_weight
  cur_network$edge_weight <- ifelse(abs(cur_network$edge_weight) > color_cutoff, sign(cur_network$edge_weight) * color_cutoff, cur_network$edge_weight)
  cur_network <- subset(cur_network, abs(edge_weight) >= cutoff)

  # make an igraph network
  cur_network <- cur_network %>%
    dplyr::rename(c(source = tf, target = gene)) %>%
    mutate(Score = sign(Cor) * Gain)

  # get a list of genes in the network
  cur_genes <- unique(cur_network$target)

  # do we only want to plot TFs?
  if (TFs_only) {
    cur_network <- subset(cur_network, target %in% unique(tf_net$tf))
  }

  # make a tidygraph
  graph <- tidygraph::as_tbl_graph(cur_network) %>%
    tidygraph::activate(nodes) %>%
    mutate(degree = centrality_degree())

  # compute the degree for each TF:
  tf_degrees <- table(tf_regulons$tf)
  tmp <- tf_degrees[names(V(graph))]
  tmp <- tmp[!is.na(tmp)]
  V(graph)[names(tmp)]$degree <- as.numeric(tmp)

  # specify the selected TFs vs TFs vs genes
  V(graph)$gene_type <- ifelse(names(V(graph)) %in% unique(tf_regulons$tf), 'TF', 'Gene')
  V(graph)$gene_type <- ifelse(names(V(graph)) %in% selected_tfs, 'selected', V(graph)$gene_type)

  # set up the graph layout
  if (igraph::is_connected(graph)) {
    n_pivots <- 250
    if (length(V(graph)) < n_pivots) {
      n_pivots <- length(V(graph)) / 2
    }
    lay <- ggraph::create_layout(graph, layout = 'sparse_stress', pivots = n_pivots)
  } else {
    lay <- ggraph::create_layout(graph, layout = 'igraph', algorithm = 'nicely')
  }

  # correct for the edge case where there are Nas
  lay[is.na(lay$x), 'x'] <- 0
  lay[is.na(lay$y), 'y'] <- 0

  # add extra info
  lay$size <- ifelse(lay$name %in% unique(tf_regulons$tf), 5, 2)

  # add the depth:
  tmp <- dplyr::left_join(lay, gene_depths, by = c('name' = 'gene'))
  lay$depth <- tmp$depth
  lay$depth <- ifelse(lay$name %in% selected_tfs, 0, lay$depth)
  lay$depth <- factor(lay$depth, levels = 0:depth)

  # shape layout:
  cur_shapes <- c(18, 17, 16)
  names(cur_shapes) <- c('selected', 'TF', 'Gene')

  # set up plotting labels
  label_tfs <- subset(cur_network, depth <= label_TFs & target %in% tf_regulons$tf) %>%
    .$target %>%
    unique()
  lay$lab <- ifelse(lay$name %in% c(selected_tfs, label_tfs, label_genes), lay$name, NA)

  # node colors
  cp <- node_colors[1:(depth + 1)]
  names(cp) <- levels(lay$depth)

  # color by depth
  p <- ggraph(lay) +
    geom_edge_link(
      aes(color = edge_weight, alpha = abs(edge_weight)),
      arrow = arrow(length = unit(1, 'mm'), type = 'closed')
    ) +
    geom_node_point(data = subset(lay, (!name %in% tf_regulons$tf)), aes(color = depth, shape = gene_type, size = degree)) +
    geom_node_point(data = subset(lay, name %in% tf_regulons$tf & !(name %in% selected_tfs)), aes(fill = depth, size = degree), color = 'black', shape = 25) +
    geom_node_point(data = subset(lay, name %in% selected_tfs), aes(color = depth, shape = gene_type, size = degree)) +
    geom_node_point(data = subset(lay, name %in% selected_tfs), aes(fill = depth, size = degree), color = 'black', shape = 23) +
    scale_edge_colour_gradient2(high = high_color, mid = mid_color, low = low_color) +
    scale_colour_manual(values = cp) +
    scale_fill_manual(values = cp) +
    scale_shape_manual(values = cur_shapes)

  # add the labels
  if (!no_labels) {
    p <- p +
      geom_node_label(aes(label = lab), repel = TRUE, max.overlaps = Inf, fontface = 'italic')
  }

  # clean up the legends
  p <- p + guides(
    edge_alpha = 'none',
    size = 'none',
    shape = 'none',
    fill = 'none'
  )
  p <- p + labs(edge_colour = 'strength')

  p
}
