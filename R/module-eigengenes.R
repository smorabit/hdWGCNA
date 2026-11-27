#' ModuleEigengenes
#'
#' Computes module eigengenes for co-expression modules
#'
#' @return seurat_obj with the module eigengenes computed for the selected wgcna experiment
#'
#' @param seurat_obj A Seurat object
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are "linear", "poisson", or "negbinom"
#' @param pc_dim Which PC to use as the module eigengene? Default to 1.
#' @param assay Assay in seurat_obj to compute module eigengenes. Default is DefaultAssay(seurat_obj)
#' @param wgcna_name name of the WGCNA experiment
#'
#' @details
#' ModuleEigengenes summarizes the gene expression signatures of entire co-expression
#' modules. This is done by performing singular value decomposition (SVD) on a
#' subset of the scaled expression matrix containing only features assigned to each module.
#' The module eigengene (ME), defined as the first dimension of the SVD matrix, retains the most variation, and we use
#' this vector as a summary of gene expression for the whole module.
#'
#' The module gene expression matrix is first scaled using the Seurat ScaleData
#' function. The user can optionally adjust for covariates of interest in this step using the
#' vars.to.regress parameter. Additionally, the module eigengenes themselves can
#' be adjusted for technical biases such as sequencing batch, dataset of origin,
#' or other factors using the Harmony algorithm with the group.by.vars parameter.
#'
#' @export
ModuleEigengenes <- function(
  seurat_obj,
  group.by.vars = NULL,
  modules = NULL,
  vars.to.regress = NULL,
  scale.model.use = 'linear',
  verbose = TRUE,
  assay = NULL,
  pc_dim = 1,
  wgcna_name = NULL, ...
) {
  # set as active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # are we going to run Harmony?
  harmonized <- !is.null(group.by.vars)

  # get the assay
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }

  me_list <- list()
  harmonized_me_list <- list()

  # re-set feature loadings:
  seurat_obj <- SetMELoadings(
    seurat_obj,
    loadings = c(''),
    harmonized = FALSE,
    wgcna_name = wgcna_name
  )
  if (harmonized) {
    seurat_obj <- SetMELoadings(
      seurat_obj,
      loadings = c(''),
      harmonized = TRUE,
      wgcna_name = wgcna_name
    )
  }

  # get modules from Seurat object, else use provided modules
  if (is.null(modules)) {
    modules <- GetModules(seurat_obj, wgcna_name)
    projected <- FALSE
  } else {
    projected <- TRUE
  }

  # get list of modules:
  mods <- levels(modules$module)
  mods_loop <- mods

  # loop over modules:
  for (cur_mod in mods_loop) {
    print(cur_mod)

    # compute module eigengenes for this module
    seurat_obj <- ComputeModuleEigengene(
      seurat_obj = seurat_obj,
      cur_mod = cur_mod,
      modules = modules,
      group.by.vars = group.by.vars,
      vars.to.regress = vars.to.regress,
      scale.model.use = scale.model.use,
      verbose = verbose,
      pc_dim = pc_dim,
      assay = assay,
      wgcna_name = wgcna_name,
      ...
    )

    # add module eigengene to ongoing list
    cur_me <- seurat_obj@reductions$ME@cell.embeddings[, pc_dim]
    me_list[[cur_mod]] <- cur_me

    # run harmony
    if (harmonized) {
      # add module eigengene to ongoing list
      cur_harmonized_me <- seurat_obj@reductions$ME_harmony@cell.embeddings[, pc_dim]
      harmonized_me_list[[cur_mod]] <- cur_harmonized_me
    }
  }

  # merge module eigengene lists into a dataframe, order modules, add to Seurat obj
  me_df <- do.call(cbind, me_list)
  if (!projected) {
    me_df <- WGCNA::orderMEs(me_df)
  }
  seurat_obj <- SetMEs(seurat_obj, me_df, harmonized = FALSE, wgcna_name)

  # merge harmonized module eigengene lists into a dataframe, add to Seurat obj
  if (!is.null(group.by.vars)) {
    hme_df <- do.call(cbind, harmonized_me_list)
    if (!projected) {
      hme_df <- WGCNA::orderMEs(hme_df)
    }
    seurat_obj <- SetMEs(seurat_obj, hme_df, harmonized = TRUE, wgcna_name)
  }

  # set module factor levels based on order
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  modules$module <- factor(
    as.character(modules$module),
    levels = mods
  )
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

  # remove temp dim reductions by setting to NULL
  seurat_obj@reductions$ME <- NULL
  seurat_obj@reductions$ME_harmony <- NULL

  # return seurat object
  seurat_obj
}


#' ComputeModuleEigengene
#'
#' Internal helper function that computes module eigengene for a single module.
#'
#' @param seurat_obj A Seurat object
#' @param cur_mod name of a module found in the modules table of the active hdWGCNA experiment (see [GetModules()])
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are "linear", "poisson", or "negbinom"
#' @param pc_dim Which PC to use as the module eigengene? Default to 1.
#' @param assay Assay in seurat_obj to compute module eigengenes. Default is DefaultAssay(seurat_obj)
#' @param wgcna_name name of the WGCNA experiment
ComputeModuleEigengene <- function(
  seurat_obj,
  cur_mod,
  modules,
  group.by.vars = NULL,
  verbose = TRUE,
  vars.to.regress = NULL,
  scale.model.use = 'linear',
  pc_dim = 1,
  assay = NULL,
  wgcna_name = NULL, ...
) {
  # set as active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get the assay
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }

  # get genes in this module:
  cur_genes <- modules %>%
    subset(module == cur_mod) %>%
    .$gene_name %>%
    as.character()

  # get the expression matrix:
  if (CheckSeurat5()) {
    X_dat <- SeuratObject::LayerData(seurat_obj, layer = 'data', assay = assay)[cur_genes, ]
  } else {
    X_dat <- Seurat::GetAssayData(seurat_obj, slot = 'data', assay = assay)[cur_genes, ]
  }

  # create seurat obj with just these genes
  cur_seurat <- CreateSeuratObject(X_dat, assay = assay, meta.data = seurat_obj@meta.data)

  if (CheckSeurat5()) {
    result <- try(cur_seurat <- SetAssayData(cur_seurat, layer = 'data', new.data = X_dat, assay = assay), silent = TRUE)
    if ('Seurat' %in% class(result)) {
      cur_seurat <- result
    } else {
      cur_seurat <- SetAssayData(cur_seurat, slot = 'data', new.data = X_dat, assay = assay)
    }
  } else {
    cur_seurat <- SetAssayData(cur_seurat, slot = 'data', new.data = X_dat, assay = assay)
  }

  # scale the subsetted expression dataset:
  if (is.null(vars.to.regress)) {
    cur_seurat <- ScaleData(cur_seurat, features = rownames(cur_seurat), model.use = scale.model.use)
  } else if (all(vars.to.regress %in% colnames(seurat_obj@meta.data))) {
    cur_seurat <- ScaleData(cur_seurat, features = rownames(cur_seurat), model.use = scale.model.use, vars.to.regress = vars.to.regress)
  } else {
    stop(paste0('Some variables specified in vars.to.regress are not found in seurat_obj@meta.data'))
  }

  # compute average expression of each gene
  if (CheckSeurat5()) {
    cur_expr <- SeuratObject::GetAssayData(cur_seurat, layer = 'data')
  } else {
    cur_expr <- Seurat::GetAssayData(cur_seurat, slot = 'data')
  }
  expr <- Matrix::t(cur_expr)
  averExpr <- Matrix::rowSums(expr) / ncol(expr)

  # run PCA with Seurat function
  cur_pca <- Seurat::RunPCA(
    cur_seurat,
    features = cur_genes,
    reduction.key = paste0('pca', cur_mod),
    verbose = verbose, ...
  )@reductions$pca
  pc <- cur_pca@cell.embeddings[, pc_dim]
  pc_loadings <- cur_pca@feature.loadings[, pc_dim]

  # correlate average expression with eigengene
  pca_cor <- cor(averExpr, pc)

  # run harmony
  if (!is.null(group.by.vars)) {
    # add this PCA as its own reduction in the seurat object
    seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
      embeddings = cur_pca@cell.embeddings,
      assay = assay
    )

    # project.dim = FALSE: hdWGCNA uses the PCA loadings, not harmony's projected loadings,
    # and projecting would require scale.data for every module gene
    harmony_args <- list(...)
    if (!('project.dim' %in% names(harmony_args))) {
      harmony_args$project.dim <- FALSE
    }
    cur_harmony <- do.call(
      harmony::RunHarmony,
      c(
        list(seurat_obj, group.by.vars = group.by.vars, reduction.use = 'ME', verbose = verbose),
        harmony_args
      )
    )@reductions$harmony
    ha <- cur_harmony@cell.embeddings[, pc_dim]
    ha_loadings <- cur_pca@feature.loadings[, pc_dim]

    if (pca_cor < 0) {
      cur_harmony@cell.embeddings[, pc_dim] <- -ha
      ha_loadings <- -ha_loadings
    }

    # add harmonized PCA as its own reduction in the seurat object
    seurat_obj@reductions$ME_harmony <- Seurat::CreateDimReducObject(
      embeddings = cur_harmony@cell.embeddings,
      assay = assay
    )

    seurat_obj <- SetMELoadings(
      seurat_obj,
      loadings = ha_loadings,
      harmonized = TRUE,
      wgcna_name = wgcna_name
    )
  }

  if (pca_cor < 0) {
    cur_pca@cell.embeddings[, pc_dim] <- -pc
    pc_loadings <- -pc_loadings
  }

  # add this PCA as its own reduction in the seurat object
  seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
    embeddings = cur_pca@cell.embeddings,
    assay = assay
  )

  seurat_obj <- SetMELoadings(
    seurat_obj,
    loadings = pc_loadings,
    harmonized = FALSE,
    wgcna_name = wgcna_name
  )

  # return seurat object
  seurat_obj
}


#' ModuleExprScore
#'
#' Computes a module score for each co-expression module using Seurat AddModuleScore or UCell.
#'
#' @return seurat_obj with module scores computed for the selected wgcna experiment
#'
#' @param seurat_obj A Seurat object
#' @param n_genes the number of genes to use for each module, ranked by kME. Setting n_genes = 'all' uses all of the genes in a module
#' @param method selected method for module scoring, valid choices are "Seurat" or "UCell"
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' ModuleExprScore provides an alternative function to ModuleEigengenes for summarizing
#' the expression level of each module. The user can choose between Seurat AddModuleScore
#' or UCell using the method parameter.
#'
#' @export
ModuleExprScore <- function(
  seurat_obj,
  n_genes = 25,
  method = 'Seurat',
  wgcna_name = NULL,
  ...
) {
  # set as active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # use all genes in the module?
  if (n_genes == 'all') {
    gene_list <- lapply(mods, function(cur_mod) {
      subset(modules, module == cur_mod) %>% .$gene_name
    })
  } else {
    gene_list <- lapply(mods, function(cur_mod) {
      cur <- subset(modules, module == cur_mod)
      cur[, c('gene_name', paste0('kME_', cur_mod))] %>%
        top_n(n_genes) %>%
        .$gene_name
    })
  }
  names(gene_list) <- mods

  # compute Module Scores with Seurat or UCell:
  if (method == 'Seurat') {
    mod_scores <- Seurat::AddModuleScore(
      seurat_obj,
      features = gene_list, ...
    )@meta.data
  } else if (method == 'UCell') {
    check_suggests('UCell', feature = 'ModuleExprScore(method = "UCell")')
    mod_scores <- UCell::AddModuleScore_UCell(
      seurat_obj,
      features = gene_list, ...
    )@meta.data
  } else {
    stop('Invalid method selection. Valid choices are Seurat, UCell')
  }

  mod_scores <- mod_scores[, (ncol(mod_scores) - length(mods) + 1):ncol(mod_scores)]

  # rename module scores:
  colnames(mod_scores) <- mods

  # order cols
  col_order <- levels(modules$module)
  col_order <- col_order[col_order != 'grey']
  mod_scores <- mod_scores[, col_order]

  # add module scores to seurat object
  seurat_obj <- SetModuleScores(seurat_obj, mod_scores, wgcna_name)

  seurat_obj
}


#' AverageModuleExpr
#'
#' Computes module eigengenes for all WGCNA co-expression modules
#'
#' @param seurat_obj A Seurat object
#' @export
AvgModuleExpr <- function(seurat_obj, n_genes = 25, wgcna_name = NULL, ...) {
  # set as active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get wgcna params
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  # datExpr for full expression dataset: the assay and layer the metacells were built from
  # (objects made with Seurat v4 store metacell_slot instead of metacell_layer)
  layer <- params$metacell_layer
  if (is.null(layer)) {
    layer <- params$metacell_slot
  }
  if (is.null(layer)) {
    layer <- 'counts'
  }
  datExpr <- SeuratObject::LayerData(seurat_obj, assay = params$metacell_assay, layer = layer)

  # use all genes in the module?
  if (n_genes == 'all') {
    gene_list <- lapply(mods, function(cur_mod) {
      subset(modules, module == cur_mod) %>% .$gene_name
    })
  } else {
    gene_list <- lapply(mods, function(cur_mod) {
      cur <- subset(modules, module == cur_mod)
      cur[, c('gene_name', paste0('kME_', cur_mod))] %>%
        top_n(n_genes) %>%
        .$gene_name
    })
  }
  names(gene_list) <- mods

  # for each module, compute average expression of genes:
  avg_mods <- t(do.call(rbind, lapply(gene_list, function(cur_genes) {
    colMeans(datExpr[cur_genes, ])
  })))

  # order cols
  col_order <- levels(modules$module)
  col_order <- col_order[col_order != 'grey']
  avg_mods <- avg_mods[, col_order]

  # add avg module expression to Seurat object
  seurat_obj <- SetAvgModuleExpr(seurat_obj, avg_mods, wgcna_name)
  seurat_obj
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
ModuleFeaturePlot <- function(
  seurat_obj, module_names = NULL, wgcna_name = NULL,
  reduction = 'umap', features = 'hMEs',
  order_points = TRUE, restrict_range = TRUE, point_size = 0.5, alpha = 1,
  label_legend = FALSE, ucell = FALSE, raster = FALSE, raster_dpi = 500,
  raster_scale = 1, plot_ratio = 1, title = TRUE
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get MEs, module data from seurat object
  if (features == 'hMEs') {
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  } else if (features == 'MEs') {
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  } else if (features == 'scores') {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else if (features == 'average') {
    MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
    restrict_range <- FALSE
  } else {
    (
      stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
    )
  }

  # override restrict_range if ucell
  if (ucell) {
    restrict_range <- FALSE
  }

  # use all modules except gray if not specified by the user
  modules <- GetModules(seurat_obj, wgcna_name)
  if (is.null(module_names)) {
    module_names <- levels(modules$module)
    module_names <- module_names[module_names != 'grey']
  }

  # get  reduction from seurat obj
  check_reduction(seurat_obj, reduction)
  umap <- seurat_obj@reductions[[reduction]]@cell.embeddings[, 1:2]
  x_name <- colnames(umap)[1]
  y_name <- colnames(umap)[2]

  # merge into one df for plotting
  plot_df <- cbind(umap, MEs) %>% as.data.frame()

  plot_list <- list()
  for (cur_mod in module_names) {
    print(cur_mod)

    # get the color for this module:
    cur_color <- modules %>%
      subset(module == cur_mod) %>%
      .$color %>%
      unique()

    # reset the range of the plot:
    plot_range <- plot_df[, cur_mod] %>% range()
    if (restrict_range) {
      if (abs(plot_range[1]) > abs(plot_range[2])) {
        plot_range[1] <- -1 * plot_range[2]
      } else {
        plot_range[2] <- -1 * plot_range[1]
      }
      plot_df[, cur_mod] <- ifelse(plot_df[, cur_mod] > plot_range[2], plot_range[2], plot_df[, cur_mod])
      plot_df[, cur_mod] <- ifelse(plot_df[, cur_mod] < plot_range[1], plot_range[1], plot_df[, cur_mod])
    }

    cur_plot_df <- plot_df[, c(colnames(umap), cur_mod)]
    colnames(cur_plot_df)[3] <- 'val'

    # order points:
    if (order_points == TRUE) {
      cur_plot_df <- cur_plot_df %>% dplyr::arrange(val)
    } else if (order_points == 'shuffle') {
      cur_plot_df <- cur_plot_df[sample(nrow(cur_plot_df)), ]
    }

    # plot with ggplot
    p <- cur_plot_df %>%
      ggplot(aes(x = .data[[x_name]], y = .data[[y_name]], color = val))

    # rasterise?
    if (raster) {
      check_suggests('ggrastr', feature = 'ModuleFeaturePlot(raster = TRUE)')
      p <- p + ggrastr::rasterise(geom_point(size = point_size, alpha = alpha), dpi = raster_dpi, scale = raster_scale)
    } else {
      p <- p + geom_point(size = point_size, alpha = alpha)
    }

    # add title and theme:
    p <- p + umap_theme() + labs(color = '')

    if (title) {
      p <- p + ggtitle(cur_mod)
    }

    # aspect ratio:
    if (is.numeric(plot_ratio)) {
      p <- p + coord_fixed(ratio = plot_ratio)
    }

    # UCell?
    if (!ucell) {
      p <- p + scale_color_gradient2(
        low = 'grey75', mid = 'grey95', high = cur_color,
        breaks = plot_range,
        labels = c('-', '+'),
        guide = guide_colorbar(ticks = FALSE, barwidth = 0.5, barheight = 4)
      )
    } else {
      p <- p + scale_color_gradient(
        low = 'grey95', high = cur_color,
        breaks = plot_range,
        labels = c('0', '+'),
        guide = guide_colorbar(ticks = FALSE, barwidth = 0.5, barheight = 4)
      )
    }
    plot_list[[cur_mod]] <- p
  }

  # return plot
  if (length(plot_list) == 1) {
    p <- plot_list[[1]]
  } else {
    p <- plot_list
  }

  p
}


#' ModuleCorrelogram
#'
#' Plot Module Eigengene correlogram
#'
#' @param seurat_obj A Seurat object
#' @param features What to plot? Can select hMEs, MEs, scores, or average
#' @keywords scRNA-seq
#' @export
ModuleCorrelogram <- function(
  seurat_obj, MEs2 = NULL,
  features = 'hMEs',
  order = 'original', method = 'ellipse',
  exclude_grey = TRUE, type = 'upper',
  tl.col = 'black', tl.srt = 45,
  sig.level = c(0.0001, 0.001, 0.01, 0.05),
  #  insig='label_sig',
  pch.cex = 0.7, col = NULL, ncolors = 200,
  wgcna_name = NULL, wgcna_name2 = NULL, ...
) {
  check_suggests('corrplot', feature = 'ModuleCorrelogram')

  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get MEs, module data from seurat object
  if (features == 'hMEs') {
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  } else if (features == 'MEs') {
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  } else if (features == 'scores') {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else if (features == 'average') {
    MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
    restrict_range <- FALSE
  } else {
    (
      stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
    )
  }

  # convert to matrix
  MEs <- as.matrix(MEs)

  # exclude grey
  if (exclude_grey) {
    MEs <- MEs[, colnames(MEs) != 'grey']
  }

  # setup color scheme:
  if (is.null(col)) {
    colfunc <- grDevices::colorRampPalette(c('seagreen', 'white', 'darkorchid1'))
    col <- colfunc(ncolors)
  }

  # perform correlation
  if (is.null(MEs2)) {
    res <- Hmisc::rcorr(x = MEs)
  } else {
    # add dataset indicator to cols/rows
    d1_names <- colnames(MEs)
    d2_names <- colnames(MEs2)
    colnames(MEs) <- paste0(d1_names, '_D1')
    colnames(MEs2) <- paste0(d2_names, '_D2')

    res <- Hmisc::rcorr(x = MEs, y = as.matrix(MEs2))

    res$r <- res$r[!grepl('_D1', colnames(res$r)), grepl('_D1', colnames(res$r))]
    colnames(res$r) <- d1_names
    rownames(res$r) <- d2_names

    res$P <- res$P[!grepl('_D1', colnames(res$P)), grepl('_D1', colnames(res$P))]
    colnames(res$P) <- d1_names
    rownames(res$P) <- d2_names
  }
  res$P[is.na(res$P)] <- 0

  # plot correlogram
  corrplot::corrplot(
    res$r,
    p.mat = res$P,
    type = type, order = order,
    method = method, tl.col = tl.col,
    tl.srt = tl.srt, sig.level = sig.level,
    # insig=insig,
    pch.cex = pch.cex,
    col = col, ...
  )
}


#' ModuleRadarPlot
#'
#' Plots the expression level (module eigengene) of each co-expression module
#' for different groups as a radar plot.
#'
#' @return ggplot object containing the ModuleRadarPlot
#'
#' @param seurat_obj A Seurat object
#' @param group.by the column name of the selected comparison in the DMEs dataframe
#' @param barcodes A list of barcodes from colnames(seurat_obj) which will be used to subset the data before plotting.
#' @param combine logical indicating whether or not to combine plots using patchwork
#' @param ncol The number of columns for the combined plot if patchwork is being used
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... additional parameters for ggradar
#' @details
#' ModuleRadarPlot visualizes the expression level (module eigengene, ME) of each co-expression module for different
#' groups in a radial coordinate system. The ME for each module is averaged for each group, scaled, and then
#' is plotted radially. The resulting plots help us to interpret which cell groups (clusters, cell types, etc) are
#' expressing each module.
#'
#' @export
ModuleRadarPlot <- function(
  seurat_obj,
  group.by = NULL,
  barcodes = NULL,
  combine = TRUE,
  ncol = 4,
  features = 'hMEs',
  wgcna_name = NULL,
  fill = TRUE,
  draw.points = FALSE,
  ... # additional params for ggradar
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  check_suggests('ggradar', feature = 'ModuleRadarPlot')

  # get seurat metadata
  meta <- seurat_obj@meta.data

  # use idents as the group
  if (is.null(group.by)) {
    cell_grouping <- Idents(seurat_obj)
  } else {
    cell_grouping <- seurat_obj@meta.data[, group.by]
    names(cell_grouping) <- colnames(seurat_obj)
  }

  # get the module info
  modules <- GetModules(seurat_obj, wgcna_name)
  mod_colors <- modules %>%
    dplyr::select(c(module, color)) %>%
    dplyr::distinct()
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get the MEs
  if (features == 'hMEs') {
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  } else if (features == 'MEs') {
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  } else if (features == 'scores') {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else if (features == 'average') {
    MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
    restrict_range <- FALSE
  } else {
    (
      stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
    )
  }

  MEs <- MEs[, colnames(MEs) != 'grey']

  # are we subsetting?
  if (!is.null(barcodes)) {
    if (!(all(barcodes %in% colnames(seurat_obj)))) {
      stop('Invalid selection for barcodes, some are not found in the colnames(seurat_obj)')
    }

    # subset:
    MEs <- MEs[barcodes, ]
    meta <- meta[barcodes, ]
    cell_grouping <- cell_grouping[barcodes]
  }

  MEs$cluster <- cell_grouping
  clusters <- as.character(unique(cell_grouping))

  # calculate the mean ME for each group
  plot_df <- MEs %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise_all(mean) %>%
    as.data.frame()

  # print(head(plot_df))

  # re-format the dataframe
  rownames(plot_df) <- plot_df$cluster
  plot_df <- dplyr::select(plot_df, -cluster)
  plot_df <- t(plot_df) %>% as.data.frame()
  plot_df[plot_df < 0] <- 0
  plot_df$group <- rownames(plot_df)
  plot_df <- plot_df[, c('group', clusters)]

  # set module factor levels
  plot_df$group <- factor(as.character(plot_df$group), levels = mods)
  plot_df <- plot_df %>%
    dplyr::arrange(group) %>%
    as.data.frame()
  colnames(plot_df) <- c('group', clusters)

  if (is.factor(cell_grouping)) {
    group_order <- levels(droplevels(cell_grouping))
  } else {
    group_order <- unique(cell_grouping)
  }

  # set the group factor levels
  plot_df <- plot_df[, c('group', group_order)]

  # make the radar plots for each module
  plot_list <- list()
  for (i in seq_len(nrow(plot_df))) {
    cur_mod <- as.character(plot_df[i, 'group'])
    cur_color <- subset(mod_colors, module == cur_mod) %>% .$color
    plot_list[[cur_mod]] <- ggradar::ggradar(
      plot_df[i, ],
      group.colours = cur_color,
      draw.points = draw.points,
      fill = fill, ...
    ) +
      Seurat::NoLegend() +
      ggtitle(cur_mod) +
      theme(
        plot.title = element_text(face = 'bold', hjust = 0.5)
      )
  }

  # combine plots with patchwork ?
  if (combine) {
    patch <- wrap_plots(plot_list, ncol)
    return(patch)
  } else {
    return(plot_list)
  }
}
