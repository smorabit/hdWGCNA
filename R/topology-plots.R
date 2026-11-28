#' ModuleTopologyHeatmap
#'
#' Plots a heatmap of the co-expression network topology of a given module.
#'
#' @return A ggplot object containing the ModuleTopologyHeatmap. If \code{return_genes = TRUE}, a list with elements \code{plot}, \code{genes} (the plotted genes in order), \code{plot_max} and \code{plot_min}.
#'
#' @param seurat_obj A Seurat object
#' @param mod the name of the co-expression module to plot
#' @param matrix specify which matrix to plot, use 'TOM' (topological overlap matrix) or 'Cor' (correlation matrix), or pass a custom square matrix where the rownames and colnames match the genes in this module
#' @param matrix_name name of the matrix plotted that will be used as the label in the plot legend
#' @param order_by order genes in this module by 'kME' (default) or by 'degree' (sum of all connections to this gene in the co-expression network)
#' @param high_color color used for high values in the heatmap, default is the module's unique color
#' @param low_color color used for low values in the heatmap, default is 'white'
#' @param raster logical indicating whether or not to rasterise the plot
#' @param raster_dpi dpi used for a rasterised plot
#' @param plot_max maximum value to plot on the heatmap, can pass a numeric value or a string indicating the quantile ('q99' would be the 99th percentile)
#' @param plot_min minimum value to plot on the heatmap, can pass a numeric value or a string indicating the quantile ('q1' would be the 1st percentile)
#' @param return_genes logical indicating whether or not to return
#' @param genes_order a character vector of genes to plot in this specific order, this option will override the order_by parameter
#' @param TOM_use  The name of the hdWGCNA experiment containing the TOM that will be used for plotting
#' @param datExpr_use  The name of the hdWGCNA experiment containing the datExpr that will be used for calculating the correlation matrix
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Additional arguments passed to \code{WGCNA::adjacency} when \code{matrix = 'Cor'}.
#' @details
#' ModuleTopologyHeatmap generates a triangular heatmap plot showing the network "topology" of a
#' specific co-expression module. Each cell in the heatmap represents a gene-gene pair, and the
#' the heatmap is colored by the strength of the connection between these two genes. By default
#' the genes in this heatmap are ordered in both the rows and the columns based on their importance
#' in the module, ranked either by eigengene-based connectivity (kME) or by network degree.
#'
#' @export
ModuleTopologyHeatmap <- function(
  seurat_obj,
  mod,
  matrix = 'TOM',
  matrix_name = NULL,
  order_by = 'kME',
  high_color = NULL,
  low_color = 'white',
  raster = TRUE,
  raster_dpi = 200,
  plot_max = 'q99',
  plot_min = 0,
  return_genes = FALSE,
  genes_order = NULL,
  TOM_use = NULL,
  datExpr_use = NULL,
  wgcna_name = NULL,
  ... # pass to WGCNA::adjacency
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  if (is.null(TOM_use)) {
    TOM_use <- wgcna_name
  }
  if (is.null(datExpr_use)) {
    datExpr_use <- wgcna_name
  }

  # is ggrastr available?
  raster_avail <- requireNamespace('ggrastr', quietly = TRUE)
  if (!raster_avail) {
    raster <- FALSE
  }

  # get the modules table
  modules <- GetModules(seurat_obj, wgcna_name)
  cur_genes <- subset(modules, module == mod) %>% .$gene_name
  hub_df <- GetHubGenes(seurat_obj, mods = mod, n_hubs = Inf, wgcna_name = wgcna_name)

  if (is.null(high_color)) {
    high_color <- subset(modules, module == mod) %>%
      .$color %>%
      unique()
  }

  # get the degree table, and subset for this module
  degrees <- GetDegrees(seurat_obj, wgcna_name) %>%
    subset(module == mod)

  if (!is.null(genes_order)) {
    genes_order <- genes_order[genes_order %in% cur_genes]
    cur_genes <- genes_order
  }

  # get the Matrix to plot:
  if (any(class(matrix) == 'matrix')) {
    mat <- matrix

    if (is.null(matrix_name)) {
      stop("Must provide matrix_name if you supply a matrix (Instead of 'TOM' or 'Cor') to the matrix argument.")
    }

    # check that it's a square matrix.
    if (!ncol(mat) == nrow(mat)) {
      stop('Invalid matrix. Must be a square matrix.')
    }

    # check that the rownames and colnames are all in the Seurat obj
    if (!all(rownames(mat) %in% rownames(seurat_obj) & all(colnames(mat) %in% rownames(seurat_obj)))) {
      stop('Invalid matrix. rownames and colnames must be found in the rownames(seurat_obj).')
    }
  } else if (matrix == 'TOM') {
    # Get the TOM
    mat <- GetTOM(seurat_obj, TOM_use)
    matrix_name <- matrix
  } else if (matrix == 'Cor') {
    # get the expression matrix:
    datExpr <- as.matrix(GetDatExpr(seurat_obj, datExpr_use))

    # calculate the unsigned adjacency matrix
    # this calculation is slow as heck (idk why)
    mat <- WGCNA::adjacency(
      datExpr,
      power = 1,
      ...
    )
    matrix_name <- matrix
  } else {
    stop("Invalid selection for matrix. Must choose 'TOM' or 'Cor', or provide a square matrix.")
  }

  # only keep the genes that are present in the matrix:
  cur_genes <- cur_genes[cur_genes %in% colnames(mat)]

  # subset tables to match the cur genes:
  hub_df <- hub_df %>% subset(gene_name %in% cur_genes)
  degrees <- degrees %>% subset(gene_name %in% cur_genes)

  # select the method to order genes:
  if (order_by == 'degree') {
    cur_genes <- as.character(degrees$gene_name)
  } else if (order_by == 'kME') {
    cur_genes <- as.character(hub_df$gene_name)
  }

  # format the matrix for plotting
  tmp <- mat[cur_genes, cur_genes]
  tmp[lower.tri(tmp)] <- 0
  plot_df <- reshape2::melt(tmp) %>% subset(Var1 != Var2)

  print('plot_df:')
  print(head(plot_df))

  # set the order for the genes
  plot_df$Var1 <- factor(plot_df$Var1, levels = cur_genes)
  plot_df$Var2 <- factor(plot_df$Var2, levels = cur_genes)

  # maximum plot values
  if (is.character(plot_max)) {
    maxquant <- as.numeric(gsub('q', '', plot_max)) / 100
    plot_max <- as.numeric(quantile(plot_df$value, maxquant))
  }
  plot_df$value <- ifelse(plot_df$value > plot_max, plot_max, plot_df$value)

  # minimum plot values
  if (is.character(plot_min)) {
    minquant <- as.numeric(gsub('q', '', plot_min)) / 100
    plot_min <- as.numeric(quantile(plot_df$value, minquant))
  }
  plot_df$value <- ifelse(plot_df$value < plot_min, plot_min, plot_df$value)

  # assemble the ggplot
  p <- plot_df %>%
    ggplot(aes(x = Var1, y = Var2, fill = value))

  # rasterise?
  if (raster) {
    p <- p + ggrastr::rasterise(geom_tile(), dpi = raster_dpi)
  } else {
    p <- p + geom_tile()
  }

  # set color
  # should I allow for gradient2?
  p <- p +
    scale_fill_gradient(low = low_color, high = high_color, limits = c(plot_min, plot_max))

  # theme
  p <- p + theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + coord_equal() + ylab('') + xlab('')

  # default plot title:
  p <- p + labs(fill = matrix_name)

  if (return_genes) {
    outs <- list()
    outs[['plot']] <- p
    outs[['genes']] <- as.character(cur_genes)
    outs[['plot_max']] <- plot_max
    outs[['plot_min']] <- plot_min
    return(outs)
  }

  # return plot
  p
}


#' ModuleTopologyBarplot
#'
#' Plots a ranked barplot of genes in a co-expression module by intramodular connectivity
#'
#' @return ggplot object containing the ModuleTopologyBarplot
#'
#' @param seurat_obj A Seurat object
#' @param mod the name of the co-expression module to plot
#' @param features specify the features to use in the barplot, 'kME' or 'degree' or 'weighted_degree' (degree scaled to 0 or 1)
#' @param plot_color color used for the bar plot, default is the module's unique color
#' @param alpha logical indicating whether or not to add opacity to the barplot based on the strength (kME or degree)
#' @param genes_order a character vector of genes to plot in this specific order, this option will override the order_by parameter
#' @param return_genes logical indicating whether or not to return
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' ModuleTopologyBarplot generates a barplot showing the intramodular connectivity of each gene
#' in a specific co-expression module. Each bar in this plot represents a single gene, and they are
#' ranked based on the strength of their connections within that particular module. A custom gene
#' ordering can be supplied, which is helpful when comparing the module topologies side by side with
#' more than one dataset.
#'
#' @export
ModuleTopologyBarplot <- function(
  seurat_obj,
  mod,
  features = 'kME', # or 'degree', or 'weighted_degree',
  plot_color = NULL, # to use the module color
  alpha = TRUE,
  genes_order = NULL,
  return_genes = FALSE,
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get the modules table
  modules <- GetModules(seurat_obj, wgcna_name)
  cur_genes <- subset(modules, module == mod) %>% .$gene_name
  hub_df <- GetHubGenes(seurat_obj, mods = mod, n_hubs = Inf, wgcna_name = wgcna_name)

  if (is.null(plot_color)) {
    plot_color <- subset(modules, module == mod) %>%
      .$color %>%
      unique()
  }

  # get the degree table, and subset for this module
  degrees <- GetDegrees(seurat_obj, wgcna_name) %>%
    subset(module == mod)

  if (features == 'kME') {
    plot_df <- hub_df %>% dplyr::rename(value = kME)
    label <- 'kME'
    plot_limits <- c(-1, 1)
  } else if (features == 'degree') {
    plot_df <- degrees %>% dplyr::rename(value = degree)
    label <- 'Degree'
    plot_limits <- c(0, max(plot_df$value))
  } else if (features == 'weighted_degree') {
    plot_df <- degrees %>% dplyr::rename(value = weighted_degree)
    label <- 'Degree'
    plot_limits <- c(0, 1)
  } else {
    stop("Invalid selection for features. Must select 'kME' or 'degree'")
  }

  # set factor levels:
  if (!is.null(genes_order)) {
    genes_order <- genes_order[genes_order %in% plot_df$gene_name]
    plot_df <- subset(plot_df, gene_name %in% genes_order)
    plot_df$gene_name <- factor(as.character(plot_df$gene_name), levels = genes_order)
    plot_df <- plot_df %>% dplyr::arrange(gene_name)
  } else {
    plot_df$gene_name <- factor(as.character(plot_df$gene_name), levels = as.character(plot_df$gene_name))
  }

  p <- plot_df %>%
    ggplot(aes(x = gene_name, y = value))

  if (alpha) {
    p <- p +
      geom_bar(aes(alpha = value), stat = 'identity', width = 1, fill = plot_color)
  } else {
    p <- p + geom_bar(stat = 'identity', width = 1, fill = plot_color)
  }

  p <- p +
    ylab(label) +
    scale_y_continuous(expand = c(0, 0), limits = plot_limits) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank()
    )

  if (return_genes) {
    return(list(p, as.character(plot_df$gene_name)))
  }

  p
}
