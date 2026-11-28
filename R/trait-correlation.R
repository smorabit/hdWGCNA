#' Module-Trait Correlation
#'
#' Computes the correlation between WGCNA modules (eigengenes or scores) and
#' biological/clinical traits (metadata columns) for single cells.
#'
#' @param seurat_obj A Seurat object containing the hdWGCNA experiment.
#' @param traits A character vector of column names in `seurat_obj@meta.data`
#'   to correlate with each module. Traits must be numeric, integer, or factor.
#'   Character vectors should be converted to factors before running this function.
#' @param group.by A string containing the name of a column in `seurat_obj@meta.data`.
#'   If provided, correlations are computed separately for each group (e.g., per cell type).
#'   If NULL (default), all cells are correlated together (global).
#' @param features Character; which feature to use to summarize each module?
#'   - "hMEs": Harmonized Module Eigengenes
#'   - "MEs": Standard Module Eigengenes.
#'   - "scores": Module Scores (Seurat `AddModuleScore`).
#' @param cor_method Character; method used for correlation. Valid choices:
#'   "pearson", "spearman", "kendall".
#' @param subset_by A string containing the name of a column to subset the data
#'   before correlation.
#' @param subset_groups A character vector specifying which groups from `subset_by`
#'   to include in the analysis.
#' @param wgcna_name The name of the hdWGCNA experiment in `seurat_obj@misc`.
#'   Default = NULL (uses active WGCNA).
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A Seurat object with the module-trait correlation results added to
#'   `seurat_obj@misc[[wgcna_name]]$mt_cor`. The results include correlation matrices,
#'   p-values, and FDR-corrected p-values for all cells and for each group specified
#'   in `group.by`.
#'
#' @details
#' This function calculates the correlation between module eigengenes and user-specified
#' traits. If a trait is a factor, it is converted to numeric based on the order of its
#' levels. The function computes p-values and False Discovery Rates (FDR) for each correlation.
#'
#' @keywords scRNA-seq
#' @export
ModuleTraitCorrelation <- function(
  seurat_obj,
  traits,
  group.by = NULL,
  features = 'hMEs',
  cor_method = 'pearson',
  subset_by = NULL,
  subset_groups = NULL,
  wgcna_name = NULL,
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get MEs, module data from seurat object
  if (features == 'hMEs') {
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  } else if (features == 'MEs') {
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  } else if (features == 'scores') {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else {
    stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
  }

  # subset?
  if (!is.null(subset_by)) {
    seurat_full <- seurat_obj
    seurat_obj <- seurat_obj[, seurat_obj@meta.data[[subset_by]] %in% subset_groups]
  }

  # Ensure MEs only contains the cells currently in seurat_obj, in the EXACT same order
  common_cells <- rownames(seurat_obj@meta.data)

  # Check if MEs contains these cells
  if (!all(common_cells %in% rownames(MEs))) {
    stop('Some cells in the Seurat object are missing from the WGCNA MEs. Please re-run ModuleEigengenes.')
  }
  MEs <- MEs[common_cells, ]

  # check if traits are in the seurat object:
  if (sum(traits %in% colnames(seurat_obj@meta.data)) != length(traits)) {
    stop(paste('Some of the provided traits were not found in the Seurat obj:', paste(traits[!(traits %in% colnames(seurat_obj@meta.data))], collapse = ', ')))
  }

  # use idents as grouping variable if not specified
  if (is.null(group.by)) {
    group.by <- 'temp_ident'
    seurat_obj$temp_ident <- Idents(seurat_obj)
  }

  # check the class of each trait provided:
  valid_types <- c('numeric', 'factor', 'integer')
  data_types <- sapply(traits, function(x) {
    class(seurat_obj@meta.data[, x])
  })

  if (!all(data_types %in% valid_types)) {
    incorrect <- traits[!(data_types %in% valid_types)]
    stop(paste0('Invalid data types for ', paste(incorrect, collapse = ', '), '. Accepted data types are numeric, factor, integer.'))
  }

  # print warnings about factor levels:
  if (any(data_types == 'factor')) {
    factor_traits <- traits[data_types == 'factor']
    for (tr in factor_traits) {
      warning(paste0('Trait ', tr, ' is a factor with levels ', paste0(levels(seurat_obj@meta.data[, tr]), collapse = ', '), '. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?'))
    }
  }

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get trait table (a data.frame even with a single trait):
  trait_df <- seurat_obj@meta.data[, traits, drop = FALSE]

  # convert factors to numeric
  if (any(data_types == 'factor')) {
    factor_traits <- traits[data_types == 'factor']
    for (tr in factor_traits) {
      trait_df[, tr] <- as.numeric(trait_df[, tr])
    }
  }

  # correlate all cells:
  cor_list <- list()
  pval_list <- list()
  fdr_list <- list()

  # correlation:
  temp <- Hmisc::rcorr(as.matrix(trait_df), as.matrix(MEs), type = cor_method)

  # get the coefficient & p-val
  cur_cor <- temp$r[traits, mods, drop = FALSE]
  cur_p <- temp$P[traits, mods, drop = FALSE]

  # compute FDR:
  p_df <- cur_p %>%
    reshape2::melt()

  if (length(traits) == 1) {
    tmp <- rep(mods, length(traits))
    tmp <- factor(tmp, levels = mods)
    tmp <- tmp[order(tmp)]

    p_df$Var1 <- traits
    p_df$Var2 <- tmp
    rownames(p_df) <- seq_len(nrow(p_df))
    p_df <- dplyr::select(p_df, c(Var1, Var2, value))
  }

  p_df <- p_df %>%
    dplyr::mutate(fdr = p.adjust(value, method = 'fdr')) %>%
    dplyr::select(c(Var1, Var2, fdr))

  # reshape to match cor & pval
  cur_fdr <- reshape2::dcast(p_df, Var1 ~ Var2, value.var = 'fdr')
  rownames(cur_fdr) <- cur_fdr$Var1
  cur_fdr <- cur_fdr[, -1]

  # add to list
  cor_list[['all_cells']] <- cur_cor
  pval_list[['all_cells']] <- cur_p
  fdr_list[['all_cells']] <- cur_fdr

  trait_df <- cbind(trait_df, seurat_obj@meta.data[, group.by])
  colnames(trait_df)[ncol(trait_df)] <- 'group'

  MEs <- cbind(as.data.frame(MEs), seurat_obj@meta.data[, group.by])
  colnames(MEs)[ncol(MEs)] <- 'group'

  if (is.factor(seurat_obj@meta.data[, group.by])) {
    group_names <- levels(seurat_obj@meta.data[, group.by])
  } else {
    group_names <- levels(as.factor(seurat_obj@meta.data[, group.by]))
  }

  # do the correlation for each group:
  combined_data <- cbind(trait_df, MEs)
  combined_data$group_label_column <- seurat_obj@meta.data[, group.by]
  split_data <- split(combined_data, combined_data$group_label_column)

  for (cur_group in names(split_data)) {
    dat <- split_data[[cur_group]]

    # Separate traits and MEs back out
    cur_traits <- dat[, traits, drop = FALSE]
    cur_MEs <- dat[, mods, drop = FALSE]

    if (nrow(dat) < 5) next

    # testing other correlation function:
    temp <- Hmisc::rcorr(as.matrix(cur_traits), as.matrix(cur_MEs), type = cor_method)
    cur_cor <- temp$r[traits, mods]
    cur_p <- temp$P[traits, mods]

    # compute FDR:
    p_df <- cur_p %>%
      reshape2::melt()

    if (length(traits) == 1) {
      tmp <- rep(mods, length(traits))
      tmp <- factor(tmp, levels = mods)
      tmp <- tmp[order(tmp)]

      p_df$Var1 <- traits
      p_df$Var2 <- tmp
      rownames(p_df) <- seq_len(nrow(p_df))
      p_df <- dplyr::select(p_df, c(Var1, Var2, value))
    }

    p_df <- p_df %>%
      dplyr::mutate(fdr = p.adjust(value, method = 'fdr')) %>%
      dplyr::select(c(Var1, Var2, fdr))

    # reshape to match cor & pval
    cur_fdr <- reshape2::dcast(p_df, Var1 ~ Var2, value.var = 'fdr')
    rownames(cur_fdr) <- cur_fdr$Var1
    cur_fdr <- cur_fdr[, -1]

    # add to list
    cor_list[[cur_group]] <- cur_cor
    pval_list[[cur_group]] <- cur_p
    fdr_list[[cur_group]] <- as.matrix(cur_fdr)
  }

  # add Module-trait correlations to the seruat object:
  mt_cor <- list(
    'cor' = cor_list,
    'pval' = pval_list,
    'fdr' = fdr_list
  )

  if (!is.null(subset_by)) {
    seurat_full <- SetModuleTraitCorrelation(seurat_full, mt_cor, wgcna_name)
    seurat_obj <- seurat_full
  } else {
    seurat_obj <- SetModuleTraitCorrelation(seurat_obj, mt_cor, wgcna_name)
  }

  seurat_obj
}


#' PlotModuleTraitCorrelation
#'
#' Plots the results of ModuleTraitCorrelation as a heatmap.
#'
#' @param seurat_obj A Seurat object
#' @param high_color color for positive correlation
#' @param mid_color color for zero correlation
#' @param low_color color for negative correlation
#' @param label which statistic to label in each cell of the heatmap: 'pval', 'fdr', or NULL (default) for no labels
#' @param label_symbol show the labels as 'stars' or as 'numeric'
#' @param plot_max maximum value of correlation to show on the colorbar
#' @param text_size size of the labels
#' @param text_color color of the text labels
#' @param text_digits how many digits to show in the text labels
#' @param combine logical determining whether to plot as one combined plot (TRUE) or to return individual plots as a list (FALSE)
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return If \code{combine = TRUE}, a patchwork object stacking one heatmap per correlation group above a module color bar; otherwise a named list of ggplot objects.
#' @keywords scRNA-seq
#' @export
PlotModuleTraitCorrelation <- function(
  seurat_obj,
  high_color = 'red',
  mid_color = 'grey90',
  low_color = 'blue',
  label = NULL,
  label_symbol = 'stars',
  plot_max = NULL,
  text_size = 2,
  text_color = 'black',
  text_digits = 3,
  combine = TRUE,
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get the module trait correlation results:
  temp <- GetModuleTraitCorrelation(seurat_obj, wgcna_name)
  cor_list <- temp$cor
  pval_list <- temp$pval
  fdr_list <- temp$fdr

  if (is.null(dim(cor_list[[1]]))) {
    stop('ModuleTraitCorrelation was run only for one trait. Heatmaps are not suggested for visualizing only one variable!')
  }

  # get module colors:
  modules <- GetModules(seurat_obj, wgcna_name)
  module_colors <- modules %>%
    dplyr::select(c(module, color)) %>%
    dplyr::distinct() %>%
    subset(module != 'grey') %>%
    dplyr::arrange(module)
  mod_colors <- module_colors$color

  # dummy variable
  module_colors$var <- 1

  # make the colorbar as its own heatmap
  module_colorbar <- module_colors %>%
    ggplot(aes(x = module, y = var, fill = module)) +
    geom_tile() +
    scale_fill_manual(values = mod_colors) +
    NoLegend() +
    RotatedAxis() +
    theme(
      plot.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )


  plot_list <- list()
  for (i in names(cor_list)) {
    cor_mat <- as.matrix(cor_list[[i]])
    pval_mat <- as.matrix(pval_list[[i]])
    fdr_mat <- as.matrix(fdr_list[[i]])

    plot_df <- reshape2::melt(cor_mat)
    colnames(plot_df) <- c('Trait', 'Module', 'cor')

    # p_df <- reshape2::melt(pval_mat)
    if (!is.null(label)) {
      if (label == 'fdr') {
        p_df <- reshape2::melt(fdr_mat)
      } else if (label == 'pval') {
        p_df <- reshape2::melt(pval_mat)
      }
      colnames(p_df) <- c('Trait', 'Module', 'pval')

      # add pval to plot_df
      plot_df$pval <- p_df$pval

      if (label_symbol == 'stars') {
        plot_df$significance <- gtools::stars.pval(plot_df$pval)
      } else if (label_symbol == 'numeric') {
        plot_df$significance <- ifelse(
          plot_df$pval <= 0.05,
          formatC(plot_df$pval, digits = text_digits), ''
        )
      } else {
        stop('Invalid input for label_symbol. Valid choices are stars or numeric.')
      }
    }

    # get limits for plot:
    if (is.null(plot_max)) {
      max_plot <- max(abs(range(plot_df$cor)))
    } else {
      max_plot <- plot_max

      # fix values outside of the specified range:
      plot_df$cor <- ifelse(abs(plot_df$cor) >= plot_max, plot_max * sign(plot_df$cor), plot_df$cor)
    }

    p <- ggplot(plot_df, aes(x = Module, y = as.numeric(Trait), fill = cor)) +
      geom_tile() +
      scale_fill_gradient2(
        limits = c(-1 * max_plot, max_plot),
        high = high_color,
        mid = mid_color,
        low = low_color,
        guide = guide_colorbar(ticks = FALSE, barwidth = 16, barheight = 0.5)
      ) +
      scale_y_continuous(
        breaks = seq_along(levels(plot_df$Trait)),
        labels = levels(plot_df$Trait),
        sec.axis = sec_axis(
          ~.,
          breaks = seq_along(levels(plot_df$Trait)),
          labels = levels(plot_df$Trait)
        )
      ) +
      RotatedAxis() +
      ylab('') +
      xlab('') +
      ggtitle(i) +
      # labs(fill = 'Correlation') +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.title = element_blank(),
        legend.position = 'bottom'
      )

    if (!is.null(label)) {
      p <- p + geom_text(label = plot_df$significance, color = text_color, size = text_size)
    }

    plot_list[[i]] <- p
  }

  if (combine) {
    # plot_list <- c(plot_list, cbar)
    # names(plot_list)[length(plot_list)] <- 'module'

    for (i in seq_along(plot_list)) {
      plot_list[[i]] <- plot_list[[i]] +
        ylab(names(plot_list)[i]) +
        theme(
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          legend.position = 'bottom',
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5)
        )
    }

    # assemble with patchwork:
    plot_list[['module']] <- module_colorbar

    out <- wrap_plots(plot_list, ncol = 1) +
      plot_layout(
        guides = 'collect',
        heights = c(rep(1, length(plot_list) - 1), 0.15)
      ) +
      plot_annotation(
        theme = theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.justification = 0.5
        )
      )

    return(out)
  } else {
    return(plot_list)
  }
}
