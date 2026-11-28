#' FindDMEs
#'
#' Function to compare expression levels of co-expression modules between two sets of cell barcodes.
#'
#' @param seurat_obj A Seurat object
#' @param barcodes1 character vector containing cell barcodes for the first group to test. Positive fold-change means up-regulated in this group.
#' @param barcodes2 character vector containing cell barcodes for the second group to test. Negative fold-change means up-regulated in this group.
#' @param features indicate whether to use "MEs" or "ModuleScores" for the comparison
#' @param harmonized logical determining whether or not to use harmonized MEs
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param add_missing logical determining whether or not to add missing modules back into the resulting dataframe with NA values.
#' @param test.use Statistical test passed to \code{Seurat::FindMarkers}. Default = 'wilcox'.
#' @param only.pos logical passed to \code{Seurat::FindMarkers}; if TRUE, only return positive results.
#' @param logfc.threshold Minimum log fold-change passed to \code{Seurat::FindMarkers}. Default = 0.
#' @param min.pct Minimum fraction of cells expressing the feature passed to \code{Seurat::FindMarkers}. Default = 0.
#' @param verbose logical passed to \code{Seurat::FindMarkers} controlling progress output.
#' @param pseudocount.use Pseudocount added before computing fold-changes, passed to \code{Seurat::FindMarkers}. Default = 0.
#' @param ... Additional parameters for the FindMarkers function. \code{latent.vars} may be given as column names of \code{seurat_obj@meta.data}.
#' @keywords scRNA-seq
#' @export
#' @return A data.frame of \code{Seurat::FindMarkers} results with one row per module and a \code{module} column.
FindDMEs <- function(
  seurat_obj,
  barcodes1,
  barcodes2,
  features = 'MEs',
  harmonized = TRUE,
  wgcna_name = NULL,
  add_missing = FALSE,
  test.use = 'wilcox',
  only.pos = FALSE,
  logfc.threshold = 0,
  min.pct = 0,
  verbose = FALSE,
  pseudocount.use = 0,
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # ensure that selected barcodes are in the seurat obj
  if (!(all(barcodes1 %in% colnames(seurat_obj)))) {
    stop('Some barcodes in barcodes1 not found in colnames(seurat_obj).')
  }
  if (!(all(barcodes2 %in% colnames(seurat_obj)))) {
    stop('Some barcodes in barcodes2 not found in colnames(seurat_obj).')
  }

  # check for overlap in the two groups of bacrodes:
  if (length(intersect(barcodes1, barcodes2)) > 0) {
    stop('Some barcodes overlap in barcodes1 and barcodes2')
  }

  # what features are we using?
  if (features == 'MEs') {
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  } else if (features == 'ModuleScores') {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else {
    stop('Invalid selection for features. Valid choices are MEs or ModuleScores.')
  }

  # remove grey module
  MEs <- MEs[, colnames(MEs) != 'grey']

  # set all negative values to zero
  MEs[MEs < 0] <- 0

  # transpose MEs so cells are columns and rows are MEs
  MEs <- t(MEs)

  # create new assay for MEs:
  ME_assay <- Seurat::CreateAssayObject(MEs)

  # run differential test on the ME assay
  DMEs <- do.call(FindMarkers, c(
    list(
      ME_assay,
      cells.1 = barcodes1,
      cells.2 = barcodes2,
      slot = 'counts',
      test.use = test.use,
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      min.pct = min.pct,
      verbose = verbose,
      pseudocount.use = pseudocount.use
    ),
    resolve_latent_vars(list(...), seurat_obj, c(barcodes1, barcodes2))
  ))
  DMEs$module <- rownames(DMEs)

  # add missing modules if specified
  if (add_missing) {
    missing_mods <- rownames(MEs)[!(rownames(MEs) %in% DMEs$module)]
    for (cur_mod in missing_mods) {
      DMEs[cur_mod, ] <- NA
      DMEs[cur_mod, 'module'] <- cur_mod
    }
  }

  DMEs
}


#' FindAllDMEs
#'
#' Function to compare expression levels of co-expression modules between two sets of cell barcodes.
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing cell grouping information
#' @param features indicate whether to use "MEs" or "ModuleScores" for the comparison
#' @param harmonized logical determining whether or not to use harmonized MEs
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param add_missing logical determining whether or not to add missing modules back into the resulting dataframe with NA values.
#' @param test.use Statistical test passed to \code{Seurat::FindMarkers}. Default = 'wilcox'.
#' @param only.pos logical passed to \code{Seurat::FindMarkers}; if TRUE, only return positive results.
#' @param logfc.threshold Minimum log fold-change passed to \code{Seurat::FindMarkers}. Default = 0.
#' @param min.pct Minimum fraction of cells expressing the feature passed to \code{Seurat::FindMarkers}. Default = 0.
#' @param verbose logical passed to \code{Seurat::FindMarkers} controlling progress output.
#' @param pseudocount.use Pseudocount added before computing fold-changes, passed to \code{Seurat::FindMarkers}. Default = 0.
#' @param ... Additional parameters for the FindMarkers function. \code{latent.vars} may be given as column names of \code{seurat_obj@meta.data}.
#' @keywords scRNA-seq
#' @export
#' @return A data.frame of \code{Seurat::FindMarkers} results for each group in \code{group.by} versus all other cells, with \code{module} and \code{group} columns.
FindAllDMEs <- function(
  seurat_obj,
  group.by,
  features = 'MEs',
  harmonized = TRUE,
  add_missing = FALSE,
  wgcna_name = NULL,
  test.use = 'wilcox',
  only.pos = FALSE,
  logfc.threshold = 0,
  min.pct = 0,
  verbose = FALSE,
  pseudocount.use = 0,
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get list of groups
  groups <- as.character(unique(seurat_obj@meta.data[[group.by]]))

  # what features are we using?
  if (features == 'MEs') {
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  } else if (features == 'ModuleScores') {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else {
    stop('Invalid selection for features. Valid choices are MEs or ModuleScores.')
  }

  # remove grey module
  MEs <- MEs[, colnames(MEs) != 'grey']

  # set all negative values to zero
  MEs[MEs < 0] <- 0

  # transpose MEs so cells are columns and rows are MEs
  MEs <- t(MEs)

  # create new assay for MEs:
  ME_assay <- Seurat::CreateAssayObject(MEs)

  DMEs_list <- list()
  for (cur_group in groups) {
    print(cur_group)

    barcodes1 <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] == cur_group]
    barcodes2 <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] != cur_group]

    # run differential test on the ME assay
    DMEs <- do.call(FindMarkers, c(
      list(
        ME_assay,
        cells.1 = barcodes1,
        cells.2 = barcodes2,
        slot = 'counts',
        test.use = test.use,
        only.pos = only.pos,
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        verbose = verbose,
        pseudocount.use = pseudocount.use
      ),
      resolve_latent_vars(list(...), seurat_obj, c(barcodes1, barcodes2))
    ))
    DMEs$module <- rownames(DMEs)
    DMEs$group <- cur_group

    # add missing modules if specified
    if (add_missing) {
      missing_mods <- rownames(MEs)[!(rownames(MEs) %in% DMEs$module)]
      for (cur_mod in missing_mods) {
        DMEs[cur_mod, ] <- NA
        DMEs[cur_mod, 'module'] <- cur_mod
      }
    }

    # add to ongoing list:
    DMEs_list[[cur_group]] <- DMEs
  }

  DMEs <- do.call(rbind, DMEs_list)
  DMEs
}


# FindMarkers() on an Assay needs latent.vars as a data.frame of per-cell values, but users
# pass metadata column names (as for a Seurat object). Convert names to that data.frame.
resolve_latent_vars <- function(dots, seurat_obj, cells) {
  lv <- dots$latent.vars
  if (!is.null(lv) && is.character(lv)) {
    missing_vars <- setdiff(lv, colnames(seurat_obj@meta.data))
    if (length(missing_vars) > 0) {
      stop(paste0(
        'latent.vars not found in seurat_obj@meta.data: ', paste0(missing_vars, collapse = ', ')
      ))
    }
    dots$latent.vars <- seurat_obj@meta.data[cells, lv, drop = FALSE]
  }
  dots
}


#' PlotDMEsVolcano
#'
#' Plotting function for the results of FindDMEs and FindAllDMEs.
#'
#' This function generates a volcano plot for differential module expression (DME) analysis results.
#' It can handle both two-group comparisons (using the output of `FindDMEs`) and one-vs-all comparisons
#' (using the output of `FindAllDMEs`).
#'
#' @param seurat_obj A Seurat object containing the WGCNA analysis in the @misc slot.
#' @param DMEs A dataframe output from FindDMEs or FindAllDMEs containing DME results.
#' @param plot_labels Logical, determines whether to plot the module labels on the volcano plot. Default is TRUE.
#' @param label_size Numeric, the size of the module labels on the plot. Default is 4.
#' @param mod_point_size Numeric, the size of the points on the volcano plot. Default is 4.
#' @param show_cutoff Logical, determines whether to plot the significance cutoff. Set this to FALSE if using facet_wrap. Default is TRUE.
#' @param wgcna_name Character, the name of the hdWGCNA experiment in the `seurat_obj@misc` slot. Default is NULL, in which case it pulls the active WGCNA experiment from `seurat_obj@misc$active_wgcna`.
#' @param xlim_range A numeric vector of length 2 specifying the x-axis limits for the log2 fold change. Default is NULL, which automatically calculates limits based on the data.
#' @param ylim_range A numeric vector of length 2 specifying the y-axis limits for the -log10(p-value). Default is NULL, which automatically calculates limits based on the data.
#' @keywords scRNA-seq, volcano plot, differential expression, WGCNA
#' @export
#' @return A ggplot object containing the volcano plot for the DME results.
#' @examples
#' # Example usage:
#' # Assuming `seurat_obj` is your Seurat object and `DMEs` is the output from FindDMEs
#' \dontrun{
#' PlotDMEsVolcano(seurat_obj, DMEs, wgcna_name = 'MG')
#' }
PlotDMEsVolcano <- function(
  seurat_obj,
  DMEs,
  plot_labels = TRUE,
  mod_point_size = 4,
  label_size = 4,
  show_cutoff = TRUE,
  wgcna_name = NULL,
  xlim_range = NULL, # New parameter to control x-axis limits
  ylim_range = NULL # New parameter to control y-axis limits
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # remove NAs:
  DMEs <- na.omit(DMEs)

  # lowest non-zero value
  lowest <- DMEs %>%
    subset(p_val_adj != 0) %>%
    top_n(-1, wt = p_val_adj) %>%
    .$p_val_adj
  DMEs$p_val_adj <- ifelse(DMEs$p_val_adj == 0, lowest, DMEs$p_val_adj)

  # fix infinite fold change
  max_fc <- max(abs(DMEs$avg_log2FC))
  max_fc <- DMEs %>%
    subset(abs(avg_log2FC) != Inf) %>%
    .$avg_log2FC %>%
    max()
  DMEs$avg_log2FC <- ifelse(DMEs$avg_log2FC == -Inf, -1 * max_fc, DMEs$avg_log2FC)
  DMEs$avg_log2FC <- ifelse(DMEs$avg_log2FC == Inf, max_fc, DMEs$avg_log2FC)

  # get modules and module colors
  modules <- GetModules(seurat_obj, wgcna_name) %>%
    subset(module != 'grey') %>%
    mutate(module = droplevels(module))
  module_colors <- modules %>%
    dplyr::select(c(module, color)) %>%
    dplyr::distinct()

  # module names
  mods <- levels(modules$module)
  mods <- mods[mods %in% DMEs$module]
  if (length(mods) == 0) {
    mods <- levels(modules$module)
    mods <- gsub('_', '-', mods) # (Probably fine to only add this line)
    mods <- mods[mods %in% DMEs$module]
  }
  mod_colors <- module_colors$color
  names(mod_colors) <- as.character(module_colors$module)

  # annotate modules with significant DME
  DMEs$anno <- ifelse(DMEs$p_val_adj < 0.05, DMEs$module, '')

  # set default x-axis limit if not provided
  if (is.null(xlim_range)) {
    xmax <- max_fc
    xlim_range <- c((-1 * xmax) - 0.1, xmax + 0.1) # Adjust as needed
  }

  # set default y-axis limit if not provided
  if (is.null(ylim_range)) {
    ymax <- max(-log10(DMEs$p_val_adj))
    ylim_range <- c(0, ymax + 1) # Adjust as needed
  }

  # plot basics
  p <- DMEs %>%
    ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), fill = module, color = module))


  if (show_cutoff) {
    p <- p +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey75', alpha = 0.8) +
      geom_rect(
        data = DMEs[1, ],
        aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)), fill = 'grey75', alpha = 0.8, color = NA
      )
  }

  # add points:
  p <- p + geom_point(size = mod_point_size, pch = 21, color = 'black')

  # label points?
  if (plot_labels) {
    p <- p + ggrepel::geom_text_repel(aes(label = anno), color = 'black', min.segment.length = 0, max.overlaps = Inf, size = label_size)
  }

  p <- p +
    scale_fill_manual(values = mod_colors) +
    scale_color_manual(values = mod_colors) +
    xlim(xlim_range) + # Apply new xlim range
    ylim(ylim_range) + # Apply new ylim range
    xlab(bquote('Average log'[2] ~ '(Fold Change)')) +
    ylab(bquote('-log'[10] ~ '(Adj. P-value)')) +
    theme(
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
      panel.grid.major = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = 'bottom'
    ) + NoLegend()

  return(p)
}


#' PlotDMEsLollipop
#'
#' Plotting function for the results of FindDMEs and FindAllDMEs
#'
#' @param seurat_obj A Seurat object
#' @param DMEs dataframe output from FindDMEs or FindAllDMEs
#' @param group.by the column name of the selected comparison in the DMEs dataframe
#' @param comparison character vector or a list of character vectors containing the comparison in the DMEs dataframe
#' @param pvalue p_value or fdr used for the comparison in the DMEs dataframe
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param avg_log2FC Name of the fold-change column in the DMEs dataframe. Currently ignored; the 'avg_log2FC' column is always used.
#' @keywords scRNA-seq
#' @export
#' @return A ggplot object if \code{group.by} is not supplied, otherwise a named list of ggplot objects, one per comparison.
#' PlotDMEsLollipop
PlotDMEsLollipop <- function(
  seurat_obj,
  DMEs,
  wgcna_name = NULL,
  group.by = NULL,
  comparison = NULL,
  pvalue = 'p_val_adj',
  avg_log2FC = 'avg_log2FC'
) {
  if (!(pvalue %in% colnames(DMEs))) {
    stop('Selected pvalue is not found in DMEs dataframe column names.')
  }

  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  modules <- GetModules(seurat_obj, wgcna_name) %>%
    subset(module != 'grey') %>%
    mutate(module = droplevels(module))

  if (!missing(group.by) && !missing(comparison)) {
    comparisons <- comparison
    if (!(all(comparisons %in% DMEs[[group.by]]))) {
      stop('Not all selected comparisons are not found in DMEs[[group.by]] or the comparison column, DMEs[[group.by]], is not correctly supplied.')
    }

    # comparisons <- unique(DMEs$comparison)
    plot_list <- list()

    for (cur_comp in comparisons) {
      print(cur_comp)

      # cur_DMEs <- DMEs %>% filter(DMEs[[group.by]] == cur_comp)
      cur_DMEs <- subset(DMEs, DMEs[[group.by]] == cur_comp)
      # provide title
      cur_title <- cur_comp
      #
      p <- PlotLollipop(modules, cur_DMEs, pvalue, avg_log2FC = 'avg_log2FC')

      p <- p + ggtitle(cur_title) + NoLegend() + geom_stripes(aes(y = module), inherit.aes = FALSE, data = cur_DMEs)

      plot_list[[cur_comp]] <- p
    }
  } else if (missing(group.by) && !missing(comparison)) { # Using && instead of & ensures that the second condition is only evaluated if the first condition is TRUE.

    stop('The group.by column is not provided in the DMEs data, and comparison cannot be found.')
  } else if (!missing(group.by) && missing(comparison)) { # Using && instead of & ensures that the second condition is only evaluated if the first condition is TRUE.

    if (!(group.by %in% names(DMEs))) {
      stop('The group.by column is not found in the DMEs data.')
    }

    comparisons <- unique(DMEs[[group.by]])

    plot_list <- list()

    for (cur_comp in comparisons) {
      print(cur_comp)

      # cur_DMEs <- DMEs %>% filter(DMEs[[group.by]] == cur_comp)
      cur_DMEs <- subset(DMEs, DMEs[[group.by]] == cur_comp)
      # provide title
      cur_title <- cur_comp

      # set plotting attributes for shape
      p <- PlotLollipop(modules, cur_DMEs, pvalue, avg_log2FC = 'avg_log2FC')

      p <- p + ggtitle(cur_title) + NoLegend() + geom_stripes(aes(y = module), inherit.aes = FALSE, data = cur_DMEs)

      plot_list[[cur_comp]] <- p
    }
  } else {
    # this is for the condition: missing(group.by) && missing(comparison)
    print('Please be aware comparison group/groups are not provided, which may casue an ERROR. PlotDMEsLollipop function will automatically assume all values are within the same group.')

    plot_list <- list()

    cur_DMEs <- DMEs

    # set plotting attributes for shape
    p <- PlotLollipop(modules, cur_DMEs, pvalue, avg_log2FC = 'avg_log2FC')

    p <- p + NoLegend() + geom_stripes(aes(y = module), inherit.aes = FALSE, data = cur_DMEs)

    plot_list <- p
  }

  return(plot_list)
}


#' PlotLollipop
#'
#' An Internal function for PlotDMEsLollipop to Plotting function for the results of FindDMEs and FindAllDMEs
#'
#' @param modules The module assignment table from \code{GetModules}, with the grey module removed.
#' @param cur_DMEs dataframe output from FindDMEs or FindAllDMEs for a single comparison
#' @param pvalue p_value or fdr used for the comparison in the DMEs dataframe
#' @param avg_log2FC Name of the fold-change column in the DMEs dataframe. Currently ignored; the 'avg_log2FC' column is always used.
#' @keywords scRNA-seq
#' @return A ggplot object
#' @examples
#' \dontrun{
#' p <- PlotLollipop(modules, DMEs, pvalue = 'p_val_adj', avg_log2FC = 'avg_log2FC')
#' }
PlotLollipop <- function(
  modules,
  cur_DMEs,
  pvalue,
  avg_log2FC = 'avg_log2FC'
) {
  #
  # cur_DMEs <- DMEs

  # set plotting attributes for shape
  cur_DMEs$shape <- ifelse(cur_DMEs[[pvalue]] < 0.05, 21, 4) # 21 cicle; 4 X
  cur_DMEs <- cur_DMEs[order(cur_DMEs[[avg_log2FC]]), ]
  cur_DMEs$module <- factor(as.character(cur_DMEs$module), levels = as.character(cur_DMEs$module))

  # add number of genes per module
  n_genes <- table(modules$module)
  cur_DMEs$n_genes <- as.numeric(n_genes[as.character(cur_DMEs$module)])

  mod_colors <- dplyr::select(modules, c(module, color)) %>% dplyr::distinct()
  cp <- mod_colors$color
  names(cp) <- mod_colors$module

  p <- cur_DMEs %>%
    ggplot(aes(y = .data$module, x = .data[[avg_log2FC]], color = .data$module)) +
    geom_vline(xintercept = 0, color = 'black') +
    geom_segment(aes(yend = .data$module, x = 0, xend = .data[[avg_log2FC]]), linewidth = 0.5, alpha = 0.3) +
    geom_point(aes(size = log(.data$n_genes))) +
    geom_point(aes(size = log(.data$n_genes)), shape = cur_DMEs$shape, color = 'black', fill = NA) +
    scale_color_manual(values = cp, guide = 'none') +
    ylab('') +
    xlab(bquote('Avg. log'[2] ~ '(Fold Change)')) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = 'plain', size = 10)
    )

  return(p)
}
