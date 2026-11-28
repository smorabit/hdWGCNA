#' ModuleEigengenes
#'
#' Computes module eigengenes for co-expression modules
#'
#' @return The input Seurat object with module eigengenes (and harmonized module eigengenes if group.by.vars is set) stored in \code{seurat_obj@misc[[wgcna_name]]}; retrieve with GetMEs().
#'
#' @param seurat_obj A Seurat object
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are "linear", "poisson", or "negbinom"
#' @param pc_dim Which PC to use as the module eigengene? Default to 1.
#' @param assay Assay in seurat_obj to compute module eigengenes. Default is DefaultAssay(seurat_obj)
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param ... Additional arguments passed to Seurat::RunPCA and, if group.by.vars is set, harmony::RunHarmony
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

  # clear loadings from any previous run; SetMELoadings() appends per module
  seurat_obj@misc[[wgcna_name]]$ME_loadings <- NULL
  if (harmonized) {
    seurat_obj@misc[[wgcna_name]]$hME_loadings <- NULL
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
    if (verbose) {
      message('Computing eigengene for module ', cur_mod)
    }

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
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param ... Additional arguments passed to Seurat::RunPCA and, if group.by.vars is set, harmony::RunHarmony
#' @return The input Seurat object with this module's eigengene stored as temporary "ME" (and "ME_harmony") reductions and its gene loadings stored in the hdWGCNA experiment.
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
#' @return The input Seurat object with module scores stored in \code{seurat_obj@misc[[wgcna_name]]}; retrieve with GetModuleScores().
#'
#' @param seurat_obj A Seurat object
#' @param n_genes the number of genes to use for each module, ranked by kME. Setting n_genes = 'all' uses all of the genes in a module
#' @param method selected method for module scoring, valid choices are "Seurat" or "UCell"
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param ... Additional arguments passed to Seurat::AddModuleScore or UCell::AddModuleScore_UCell, depending on method
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
#' Computes the average expression of the top genes in each co-expression module,
#' using the assay and layer that the metacells were constructed from.
#'
#' @param seurat_obj A Seurat object.
#' @param n_genes the number of genes to use for each module, ranked by kME. Setting n_genes = 'all' uses all of the genes in a module. Default = 25
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param ... Not currently used.
#' @return The input Seurat object with average module expression stored in \code{seurat_obj@misc[[wgcna_name]]}; retrieve with GetAvgModuleExpr().
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
#' @param seurat_obj A Seurat object.
#' @param module_names Character vector of modules to plot. Default = NULL to plot all modules except grey.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param reduction Name of the dimensionality reduction to plot. Default = 'umap'
#' @param features What to plot? Can select hMEs, MEs, scores, or average
#' @param order_points How to order points: TRUE (highest values on top), FALSE, or "shuffle". Default = TRUE
#' @param restrict_range Logical; if TRUE, clip the color scale to be symmetric around zero. Default = TRUE
#' @param point_size Size of the points. Default = 0.5
#' @param alpha Transparency of the points. Default = 1
#' @param label_legend Not currently used.
#' @param ucell Logical; set to TRUE when plotting UCell scores to use a one-sided color scale. Default = FALSE
#' @param raster Logical; if TRUE, rasterise the points with ggrastr. Default = FALSE
#' @param raster_dpi Resolution of the rasterised points. Default = 500
#' @param raster_scale Scaling factor for rasterised points. Default = 1
#' @param plot_ratio Aspect ratio passed to ggplot2::coord_fixed; use a non-numeric value to skip. Default = 1
#' @param title Logical; if TRUE, add the module name as the plot title. Default = TRUE
#' @return A ggplot object if one module is plotted, otherwise a named list of ggplot objects (one per module).
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
#' Plot the correlation between modules as a ggplot heatmap.
#'
#' Correlations and their p-values are computed with \code{Hmisc::rcorr}. Each
#' cell shows the correlation as a color (and, for \code{method = 'circle'}, as
#' the circle size), with stars marking significant correlations.
#'
#' @param seurat_obj A Seurat object.
#' @param MEs2 An optional matrix or data frame of module features (cells by modules) from a second dataset, measured in the same cells, to correlate against. Default = NULL
#' @param features What to plot? Can select hMEs, MEs, scores, or average
#' @param order Order of the modules: 'original' (module order) or 'hclust' (hierarchical clustering of the correlation matrix). Default = 'original'
#' @param method 'circle' (default; circle size and color show the correlation) or 'tile' (colored squares).
#' @param exclude_grey Logical; if TRUE, exclude the grey module. Default = TRUE
#' @param type Which part of the matrix to display: 'upper' (default), 'lower' or 'full'. Ignored when \code{MEs2} is given.
#' @param tl.col Color of the module labels. Default = 'black'
#' @param tl.srt Rotation angle of the x-axis module labels. Default = 45
#' @param sig.level Increasing p-value thresholds marked with three, two and one star. Default = c(0.001, 0.01, 0.05); NULL to hide stars.
#' @param pch.cex Size of the significance stars. Default = 3
#' @param col Vector of colors for the color scale, from -1 to 1. Default = NULL to use a seagreen-white-purple gradient.
#' @param ncolors Not currently used; kept for compatibility.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param wgcna_name2 Not currently used.
#' @param ... Ignored; previously passed to corrplot::corrplot.
#' @return A ggplot object.
#' @keywords scRNA-seq
#' @section Deprecated:
#' This function is deprecated and will be removed in the next release.
#' @keywords internal
#' @export
ModuleCorrelogram <- function(
  seurat_obj, MEs2 = NULL,
  features = 'hMEs',
  order = c('original', 'hclust'),
  method = c('circle', 'tile'),
  exclude_grey = TRUE, type = c('upper', 'lower', 'full'),
  tl.col = 'black', tl.srt = 45,
  sig.level = c(0.001, 0.01, 0.05),
  pch.cex = 3, col = NULL, ncolors = 200,
  wgcna_name = NULL, wgcna_name2 = NULL, ...
) {
  .Deprecated(msg = paste0(
    'ModuleCorrelogram() is deprecated and will be removed in the next release of hdWGCNA.'
  ))
  order <- match.arg(order)
  method <- match.arg(method)
  type <- match.arg(type)
  if (...length() > 0) {
    warning('ModuleCorrelogram no longer uses corrplot; additional arguments in ... are ignored.')
  }
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get MEs, module data from seurat object
  MEs <- switch(features,
    hMEs = GetMEs(seurat_obj, TRUE, wgcna_name),
    MEs = GetMEs(seurat_obj, FALSE, wgcna_name),
    scores = GetModuleScores(seurat_obj, wgcna_name),
    average = GetAvgModuleExpr(seurat_obj, wgcna_name),
    stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
  )
  MEs <- as.matrix(MEs)
  if (exclude_grey) {
    MEs <- MEs[, colnames(MEs) != 'grey', drop = FALSE]
  }

  # correlation and p-values
  if (is.null(MEs2)) {
    res <- Hmisc::rcorr(x = MEs)
    r <- res$r
    p <- res$P
  } else {
    MEs2 <- as.matrix(MEs2)
    if (nrow(MEs2) != nrow(MEs)) {
      stop('MEs2 must have one row per cell, like the module features of seurat_obj.')
    }
    res <- Hmisc::rcorr(x = MEs, y = MEs2)
    d1 <- seq_len(ncol(MEs))
    d2 <- ncol(MEs) + seq_len(ncol(MEs2))
    r <- res$r[d2, d1, drop = FALSE]
    p <- res$P[d2, d1, drop = FALSE]
    rownames(r) <- rownames(p) <- colnames(MEs2)
    colnames(r) <- colnames(p) <- colnames(MEs)
    type <- 'full'
  }

  # order
  row_order <- rownames(r)
  col_order <- colnames(r)
  if (order == 'hclust' && is.null(MEs2)) {
    hc <- stats::hclust(stats::as.dist(1 - r))
    row_order <- col_order <- rownames(r)[hc$order]
  }

  plot_df <- data.frame(
    row = rep(rownames(r), times = ncol(r)),
    col = rep(colnames(r), each = nrow(r)),
    r = as.vector(r),
    p = as.vector(p)
  )
  plot_df$row <- factor(plot_df$row, levels = rev(row_order))
  plot_df$col <- factor(plot_df$col, levels = col_order)

  # keep one triangle (including the diagonal)
  if (type != 'full') {
    i <- match(plot_df$row, row_order)
    j <- match(plot_df$col, col_order)
    plot_df <- plot_df[if (type == 'upper') i <= j else i >= j, ]
  }

  # significance stars
  plot_df$stars <- ''
  if (!is.null(sig.level)) {
    sig.level <- sort(sig.level)
    for (k in rev(seq_along(sig.level))) {
      hit <- !is.na(plot_df$p) & plot_df$p < sig.level[k]
      plot_df$stars[hit] <- strrep('*', length(sig.level) - k + 1)
    }
  }

  if (is.null(col)) {
    col <- c('seagreen', 'white', 'darkorchid1')
  }

  p_out <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$col, y = .data$row))
  if (method == 'circle') {
    p_out <- p_out +
      ggplot2::geom_tile(fill = NA, color = 'grey90') +
      ggplot2::geom_point(
        ggplot2::aes(size = abs(.data$r), fill = .data$r),
        shape = 21, color = 'grey40', stroke = 0.2
      ) +
      ggplot2::scale_size_area(max_size = 10, limits = c(0, 1), guide = 'none')
  } else {
    p_out <- p_out + ggplot2::geom_tile(ggplot2::aes(fill = .data$r), color = 'white')
  }
  # stars sit at the top of the cell so they do not hide small circles
  star_nudge <- if (method == 'circle') 0.32 else 0
  p_out +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$stars),
      size = pch.cex, vjust = 0.75, position = ggplot2::position_nudge(y = star_nudge)
    ) +
    ggplot2::scale_fill_gradientn(colors = col, limits = c(-1, 1), name = 'Correlation') +
    ggplot2::scale_x_discrete(position = 'top') +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = tl.col),
      axis.text.x.top = ggplot2::element_text(angle = tl.srt, hjust = 0, vjust = 0)
    )
}


#' ModuleRadarPlot
#'
#' Plots the expression level (module eigengene) of each co-expression module
#' for different groups as a radar plot.
#'
#' @return A patchwork object of radar plots if combine = TRUE, otherwise a named list of ggplot objects (one per module).
#'
#' @param seurat_obj A Seurat object
#' @param group.by Column in \code{seurat_obj@meta.data} defining the groups (spokes of the radar). Default = NULL to use \code{Idents(seurat_obj)}.
#' @param barcodes A list of barcodes from colnames(seurat_obj) which will be used to subset the data before plotting.
#' @param combine logical indicating whether or not to combine plots using patchwork
#' @param ncol The number of columns for the combined plot if patchwork is being used
#' @param features What to plot? Can select hMEs, MEs, scores, or average
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param fill Logical; if TRUE, fill the radar polygons. Default = TRUE
#' @param draw.points Logical; if TRUE, draw points at each vertex. Default = FALSE
#' @param ... Ignored; previously passed to ggradar::ggradar.
#' @details
#' ModuleRadarPlot visualizes the expression level (module eigengene, ME) of each co-expression module for different
#' groups in a radial coordinate system. The ME for each module is averaged for each group; negative averages are set
#' to zero and the averages are scaled by the module's maximum, so the group with the highest expression reaches the
#' outer ring. The resulting plots help us to interpret which cell groups (clusters, cell types, etc) are
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
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  if (...length() > 0) {
    warning('ModuleRadarPlot no longer uses ggradar; additional arguments in ... are ignored.')
  }

  # groups
  if (is.null(group.by)) {
    cell_grouping <- Idents(seurat_obj)
  } else {
    cell_grouping <- seurat_obj@meta.data[[group.by]]
  }
  names(cell_grouping) <- colnames(seurat_obj)

  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- setdiff(levels(modules$module), 'grey')
  mod_colors <- as.character(modules$color[match(mods, modules$module)])

  MEs <- switch(features,
    hMEs = GetMEs(seurat_obj, TRUE, wgcna_name),
    MEs = GetMEs(seurat_obj, FALSE, wgcna_name),
    scores = GetModuleScores(seurat_obj, wgcna_name),
    average = GetAvgModuleExpr(seurat_obj, wgcna_name),
    stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
  )
  MEs <- as.matrix(MEs)[, mods, drop = FALSE]

  if (!is.null(barcodes)) {
    if (!all(barcodes %in% colnames(seurat_obj))) {
      stop('Invalid selection for barcodes, some are not found in the colnames(seurat_obj)')
    }
    MEs <- MEs[barcodes, , drop = FALSE]
    cell_grouping <- cell_grouping[barcodes]
  }
  cell_grouping <- droplevels(factor(cell_grouping))
  groups <- levels(cell_grouping)

  # mean per group (rows) and module (columns), negatives to 0, scaled to the module max
  avg <- rowsum(MEs, cell_grouping) / as.vector(table(cell_grouping))
  avg[avg < 0] <- 0
  col_max <- apply(avg, 2, max)
  avg <- sweep(avg, 2, ifelse(col_max > 0, col_max, 1), '/')

  # radar geometry: one spoke per group, clockwise from the top
  n_groups <- length(groups)
  angle <- pi / 2 - 2 * pi * (seq_len(n_groups) - 1) / n_groups
  rings <- c(0.25, 0.5, 0.75, 1)
  ring_df <- data.frame(
    ring = rep(rings, each = n_groups + 1),
    x = rep(rings, each = n_groups + 1) * cos(c(angle, angle[1])),
    y = rep(rings, each = n_groups + 1) * sin(c(angle, angle[1]))
  )
  spoke_df <- data.frame(x = 0, y = 0, xend = cos(angle), yend = sin(angle))
  label_df <- data.frame(
    label = groups,
    x = 1.18 * cos(angle),
    y = 1.18 * sin(angle),
    hjust = ifelse(abs(cos(angle)) < 0.2, 0.5, ifelse(cos(angle) > 0, 0, 1))
  )

  plot_list <- lapply(seq_along(mods), function(i) {
    r <- avg[groups, mods[i]]
    poly_df <- data.frame(x = r * cos(angle), y = r * sin(angle))
    p <- ggplot2::ggplot() +
      ggplot2::geom_path(
        data = ring_df,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$ring),
        color = 'grey80', linewidth = 0.3
      ) +
      ggplot2::geom_segment(
        data = spoke_df,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        color = 'grey80', linewidth = 0.3
      ) +
      ggplot2::geom_polygon(
        data = poly_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        fill = if (fill) mod_colors[i] else NA,
        alpha = 0.3, color = mod_colors[i], linewidth = 0.8
      ) +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, hjust = .data$hjust),
        size = 3
      )
    if (draw.points) {
      p <- p + ggplot2::geom_point(data = poly_df, ggplot2::aes(x = .data$x, y = .data$y), color = mod_colors[i])
    }
    p +
      ggplot2::coord_fixed(clip = 'off', xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
      ggplot2::ggtitle(mods[i]) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = 'bold', hjust = 0.5),
        plot.margin = ggplot2::margin(10, 30, 10, 30)
      )
  })
  names(plot_list) <- mods

  if (combine) {
    return(patchwork::wrap_plots(plot_list, ncol = ncol))
  }
  plot_list
}
