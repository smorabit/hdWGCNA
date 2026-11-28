#' ModuleConnectivity
#'
#' Computes eigengene-based connectivity (kME) based on module eigengenes and the
#' feature expression in selected cell populations / spatial regions.
#'
#' @return The input Seurat object with kME columns added to the modules table in \code{seurat_obj@misc[[wgcna_name]]} (retrieve with GetModules() or GetHubGenes()), and intramodular degrees if the TOM is available (retrieve with GetDegrees()).
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing grouping info, ie clusters or celltypes
#' @param group_name name of the group(s) in group.by to use for kME calculation
#' @param corFnc character string specifying the function to be used to calculate co-expression similarity. Defaults to bicor. Any function returning values
#' @param corOptions 	character string specifying additional arguments to be passed to the function given by corFnc. Use "use = 'p', method = 'spearman'" to obtain Spearman correlation. Use "use = 'p'" to obtain Pearson correlation.
#' @param harmonized logical determining whether to use harmonized MEs for kME calculation
#' @param assay Assay in seurat_obj containing expression information.
#' @param slot Slot in seurat_obj, default to normalized 'data' slot. Slot is used with Seurat v4 instead of layer.
#' @param layer Layer in seurat_obj, default to normalized 'data' layer. Layer is used with Seurat v5 instead of slot.
#' @param sparse logical indicating whether or not to run the correlation using a sparse matrix.
#' @param reassign_modules logical indicating whether or not to reassign genes to different co-expression modules if their kME value in the assigned module is negative.
#' @param TOM_use name of the hdWGCNA experiment that contains the TOM which will be used to compute the intramodular degree of each gene. If the TOM is not found, this calculation will be skipped.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param ... Additional arguments passed to WGCNA::signedKME (only used when sparse = FALSE)
#' @details
#' ModuleConnectivity computes the eigengene-based connectivity (kME) of each feature
#' and each module. This is done by correlating the gene expression signal with
#' the module eigengenes for each module. Features with high kME values have greater
#' network connectivity within that module, and we can identify module hub genes
#' by taking the top genes for each module by kME values.
#'
#' We recommend using the group.by and group_name parameters that were previously
#' used in the SetDatExpr function, so only the relevant cell populations / spatial
#' regions are considered for computing kMEs.
#'
#' Depending on the parameters for network analysis, sometimes a feature has a negative kME
#' in the module that it was originally assigned. This discrepancy arises since
#' the network construction is done on the metacell/metaspot matrix but the kME
#' calculation is done on the cell/spot matrix. By default, we account for this
#' by reassigning features with negative kMEs in their assigned module to the
#' module that had the highest kME for that feature.
#'
#' @export
ModuleConnectivity <- function(
  seurat_obj,
  group.by = NULL,
  group_name = NULL,
  corFnc = 'bicor',
  corOptions = "use='p'",
  harmonized = TRUE,
  assay = NULL,
  slot = 'data',
  layer = 'data',
  sparse = TRUE,
  reassign_modules = TRUE,
  TOM_use = NULL,
  wgcna_name = NULL,
  ...
) {
  # set as active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  if (is.null(TOM_use)) {
    TOM_use <- wgcna_name
  }

  # get module df, wgcna genes, TOM, and wgcna params:
  modules <- GetModules(seurat_obj, wgcna_name)
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  genes_use <- as.character(modules$gene_name)
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  # try to load the TOM
  TOM <- tryCatch(
    {
      GetTOM(seurat_obj, TOM_use)
    },
    error = function(e) {
      warning('TOM not found')
    }
  )

  # get the assay
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }

  # select which cells to use for this analysis
  if (!is.null(group.by)) {
    cells.use <- seurat_obj@meta.data %>%
      subset(get(group.by) %in% group_name) %>%
      rownames()
    MEs <- MEs[cells.use, ]
  } else {
    cells.use <- colnames(seurat_obj)
  }

  # get expression data:
  if (CheckSeurat5()) {
    exp_mat <- SeuratObject::LayerData(
      seurat_obj,
      assay = assay, layer = layer
    )[genes_use, cells.use]
  } else {
    exp_mat <- Seurat::GetAssayData(
      seurat_obj,
      assay = assay, slot = slot
    )[genes_use, cells.use]
  }

  # are we using the sparse correlation? (faster!)
  if (sparse) {
    kMEs <- corSparse(
      X = Matrix::t(exp_mat),
      Y = as.matrix(MEs)
    )
    rownames(kMEs) <- genes_use
    kMEs <- as.data.frame(kMEs)
  } else {
    datExpr <- t(as.matrix(exp_mat))
    kMEs <- WGCNA::signedKME(
      datExpr,
      MEs,
      outputColumnName = 'kME',
      corFnc = corFnc,
      corOptions = corOptions,
      ...
    )
  }

  # add module color to the kMEs table
  modules <- modules[, 1:3]
  mods <- levels(modules$module)
  colnames(kMEs) <- colnames(MEs)
  kMEs <- kMEs[, mods]
  colnames(kMEs) <- paste0('kME_', colnames(kMEs))
  kMEs <- cbind(modules, kMEs)

  # update the modules table in the Seurat object:
  seurat_obj <- SetModules(seurat_obj, kMEs, wgcna_name)

  # reassign modules for genes with negative kME values
  if (reassign_modules && length(mods) > 2) {
    seurat_obj <- ReassignModules(seurat_obj, harmonized = harmonized, wgcna_name = wgcna_name)
  }

  # compute the intramodular degree for each gene
  if ('matrix' %in% class(TOM)) {
    degrees <- data.frame()
    for (cur_mod in mods) {
      # get genes for this module
      cur_genes <- modules %>%
        subset(module == cur_mod) %>%
        .$gene_name
      cur_genes <- cur_genes[cur_genes %in% rownames(TOM)]

      # subset TOM for genes in this module
      cur_TOM <- TOM[cur_genes, cur_genes]
      cur_TOM[lower.tri(cur_TOM)] <- 0

      # cast from from matrix to long dataframe, remove self-links
      cur_TOM_df <- reshape2::melt(cur_TOM) %>% subset(Var1 != Var2)

      # calculate degree for each gene and arrange by degree
      cur_degrees <- cur_TOM_df %>%
        dplyr::group_by(Var1) %>%
        dplyr::summarise(degree = sum(value)) %>%
        dplyr::mutate(weighted_degree = degree / max(degree)) %>%
        dplyr::arrange(-degree) %>%
        dplyr::rename(gene_name = Var1)

      cur_degrees$module <- cur_mod
      degrees <- rbind(degrees, cur_degrees)
    }
    degrees <- degrees %>% dplyr::select(c(gene_name, module, degree, weighted_degree))
    degrees <- as.data.frame(degrees)
    degrees$module <- factor(as.character(degrees$module), levels = levels(modules$module))
    degrees$gene_name <- as.character(degrees$gene_name)

    # update the degree df in the Seurat object:
    seurat_obj <- SetDegrees(seurat_obj, degrees, wgcna_name)
  }

  # return seurat obj
  seurat_obj
}


#' PlotKMEs
#'
#' Plotting function to show genes by kME value in each module
#'
#' @param seurat_obj A Seurat object
#' @param n_hubs number of hub genes to display
#' @param text_size controls the size of the hub gene text
#' @param ncol number of columns to display individual plots
#' @param plot_widths the relative width between the kME rank plot and the hub gene text
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A patchwork object with one kME rank plot and hub gene label panel per module (excluding grey).
#' @keywords scRNA-seq
#' @export
PlotKMEs <- function(
  seurat_obj,
  n_hubs = 10,
  text_size = 2,
  ncol = 5,
  plot_widths = c(3, 2),
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 'grey')
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']
  mod_colors <- modules %>%
    subset(module %in% mods) %>%
    dplyr::select(c(module, color)) %>%
    dplyr::distinct()

  # get hub genes:
  hub_df <- GetHubGenes(seurat_obj, n_hubs = n_hubs, wgcna_name = wgcna_name)

  plot_list <- lapply(mods, function(x) {
    cur_color <- subset(mod_colors, module == x) %>% .$color
    cur_df <- subset(hub_df, module == x)
    top_genes <- cur_df %>%
      dplyr::top_n(n_hubs, wt = kME) %>%
      .$gene_name
    p <- cur_df %>% ggplot(aes(x = reorder(gene_name, kME), y = kME)) +
      geom_bar(stat = 'identity', width = 1, color = cur_color, fill = cur_color) +
      ggtitle(x) +
      # xlab(paste0('kME_', x)) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.line.x = element_blank()
      )
    p_anno <- ggplot() +
      annotate(
        'label',
        x = 0,
        y = 0,
        label = paste0(top_genes, collapse = '\n'),
        size = text_size,
        fontface = 'italic',
        label.size = 0
      ) +
      theme_void()
    patch <- p + p_anno + plot_layout(widths = plot_widths)
    patch
  })

  wrap_plots(plot_list, ncol = ncol)
}


#' Plots gene expression of hub genes as a heatmap
#'
#' This function makes an expression heatmap of the top n hub genes per module
#' using Seurat's DoHeatmap, and then assembles them all into one big heatmap.
#'
#'
#'
#' @param seurat_obj A Seurat object.
#' @param n_hubs Number of hub genes per module to plot, ranked by kME. Default = 10
#' @param n_cells Maximum number of cells to sample from each group. Default = 200
#' @param group.by Column in seurat_obj@meta.data used to group cells. Default = NULL to use Idents(seurat_obj).
#' @param module_names Character vector of modules to plot. Default = NULL to plot all modules except grey.
#' @param combine Logical; if TRUE, combine the heatmaps with a module colorbar into one patchwork; if FALSE return a list of heatmaps. Default = TRUE
#' @param draw.lines Logical; whether to draw lines separating groups, passed to Seurat::DoHeatmap. Default = TRUE
#' @param disp.min Minimum display value for scaled expression. Default = -2.5
#' @param disp.max Maximum display value for scaled expression. Default = 2.5
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A patchwork object if combine = TRUE, otherwise a list of ggplot heatmaps (one per module).
#' @keywords scRNA-seq
#' @export
DoHubGeneHeatmap <- function(
  seurat_obj,
  n_hubs = 10,
  n_cells = 200,
  group.by = NULL,
  module_names = NULL,
  combine = TRUE, # returns a list of individual heatmaps if FALSE
  draw.lines = TRUE,
  disp.min = -2.5, disp.max = 2.5, # cutoff expression values
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # use idents as grouping variable if not specified
  if (is.null(group.by)) {
    group.by <- 'temp_ident'
    seurat_obj$temp_ident <- Idents(seurat_obj)
  }

  # drop if there are missing levels:
  seurat_obj@meta.data[[group.by]] <- droplevels(seurat_obj@meta.data[[group.by]])

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  modules <- modules %>%
    subset(module != 'grey') %>%
    mutate(module = droplevels(module))
  mods <- levels(modules$module)

  if (!is.null(module_names)) {
    print('here')
    mods <- module_names
    modules <- modules %>% subset(module %in% mods)
  }

  # get table of module names & colors
  mod_colors <- modules %>%
    dplyr::select(c(module, color)) %>%
    dplyr::distinct()

  # get hub genes:
  # hub_list <- lapply(mods, function(cur_mod){
  #   cur <- subset(modules, module == cur_mod)
  #   cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
  #     top_n(n_hubs) %>% .$gene_name
  # })
  hub_list <- lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur <- cur[, c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(n_hubs)
    colnames(cur)[2] <- 'var'
    cur %>%
      arrange(desc(var)) %>%
      .$gene_name
  })
  names(hub_list) <- mods
  print(hub_list)

  seurat_obj$barcode <- colnames(seurat_obj)
  temp <- table(seurat_obj@meta.data[[group.by]])

  # sample cells
  df <- data.frame()
  for (i in seq_along(temp)) {
    if (temp[[i]] < n_cells) {
      cur_df <- seurat_obj@meta.data %>% subset(get(group.by) == names(temp)[i])
    } else {
      cur_df <- seurat_obj@meta.data %>%
        subset(get(group.by) == names(temp)[i]) %>%
        sample_n(n_cells)
    }
    df <- rbind(df, cur_df)
  }

  # make sampled seurat obj for plotting:
  seurat_plot <- seurat_obj %>% subset(barcode %in% df$barcode)

  plot_list <- list()
  for (i in seq_along(hub_list)) {
    print(i)
    cur_mod <- names(hub_list)[i]
    print(i)
    print(hub_list[[i]])
    print(i)

    if (i == 1) {
      plot_list[[i]] <- DoHeatmap(
        seurat_plot,
        features = hub_list[[i]],
        group.by = group.by,
        raster = TRUE, slot = 'scale.data',
        disp.min = disp.min, disp.max = disp.max,
        label = FALSE,
        group.bar = FALSE,
        draw.lines = draw.lines
      )
    } else {
      plot_list[[i]] <- DoHeatmap(
        seurat_plot,
        features = hub_list[[i]],
        group.by = group.by,
        raster = TRUE, slot = 'scale.data',
        group.bar.height = 0,
        label = FALSE, group.bar = FALSE,
        draw.lines = draw.lines,
        disp.min = disp.min, disp.max = disp.max
      ) + NoLegend()
    }
    print(i)
    # margin:
    plot_list[[i]] <- plot_list[[i]] +
      theme(
        plot.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(face = 'italic')
      ) + scale_y_discrete(position = 'right')
    print(i)
  }

  # useful for ratios on colorbars
  n_total_cells <- ncol(seurat_plot)
  width_cbar <- n_total_cells / 50


  # module colorbar
  mod_colors$value <- n_hubs
  mod_colors$dummy <- 'colorbar'
  cbar_list <- list()
  for (i in seq_len(nrow(mod_colors))) {
    cbar_list[[i]] <- mod_colors[i, ] %>% ggplot(aes(y = value, x = dummy)) +
      geom_bar(position = 'stack', stat = 'identity', fill = mod_colors[i, ]$color) +
      umap_theme() +
      theme(
        plot.margin = margin(0, 0, 0, 0)
      )
  }
  p_cbar <- wrap_plots(cbar_list, ncol = 1)

  if (combine) {
    out <- wrap_plots(plot_list, ncol = 1) + plot_layout(guides = 'collect')
    out <- (p_cbar | out) + plot_layout(widths = c(width_cbar, n_total_cells))
  } else {
    out <- plot_list
  }

  out
}
