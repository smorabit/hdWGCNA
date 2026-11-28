#' ModulePreservation
#'
#' Computes module preservation statistics in a query dataset for a given reference dataset
#'
#' @param seurat_obj A Seurat object
#' @param seurat_ref A Seurat object serving as the reference for the module preservation analysis
#' @param name The name to give the module preservation analysis.
#' @param n_permutations Number of permutations for the module preservation test.
#' @param parallel logical determining whether to run preservation analysis in parallel
#' @param seed random seed for the permutation analysis.
#' @param gene_mapping a dataframe containing gene name mappings between the query and the reference dataset. One column should have the gene name in the query dataset, and anotehr column should have the corresponding gene name in the reference dataset.
#' @param genome_ref_col the column name containing the gene names for the reference dataset
#' @param genome_query_col the column name containing the gene names for the query dataset
#' @param return_raw if TRUE, returns the module preservation statistics, else returns seurat_obj with the stats added to the hdWGCNA experiment.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param wgcna_name_ref The name of the hdWGCNA experiment in the seurat_ref@misc slot. Defaults to the active experiment in seurat_ref.
#' @param ... Additional arguments passed to WGCNA::modulePreservation
#' @return If return_raw = TRUE, the raw output list of WGCNA::modulePreservation; otherwise seurat_obj with the observed and Z-score preservation statistics stored under name in \code{seurat_obj@misc[[wgcna_name]]} (retrieve with GetModulePreservation()).
#' @details
#' ModulePreservation performs a statistical test to assess the preservation of co-expression modules identified
#' in one dataset in an independent dataset. This method is originally described by Langfelder et al.
#' in the 2011 paper "Is My Network Module Preserved and Reproducible?". This method can be used to
#' assess biological differences in networks, as well as technical differences / reproducibility across different
#' batches.
#'
#' @export
ModulePreservation <- function(
  seurat_obj,
  seurat_ref,
  name,
  n_permutations = 500,
  parallel = FALSE,
  seed = 12345,
  gene_mapping = NULL,
  genome_ref_col = NULL,
  genome_query_col = NULL,
  return_raw = FALSE,
  wgcna_name = NULL,
  wgcna_name_ref = NULL,
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  if (is.null(wgcna_name_ref)) {
    wgcna_name_ref <- seurat_ref@misc$active_wgcna
  }
  CheckWGCNAName(seurat_ref, wgcna_name_ref)

  # get datExpr for reference and query:
  datExpr_ref <- GetDatExpr(seurat_ref, wgcna_name_ref)
  datExpr_query <- GetDatExpr(seurat_obj, wgcna_name)

  # change the gene names to match:
  if (!is.null(gene_mapping)) {
    gene_match <- match(colnames(datExpr_query), gene_mapping[, genome_query_col])
    gene_mapping <- na.omit(gene_mapping[gene_match, ])
    colnames(datExpr_query) <- gene_mapping[, genome_ref_col]
  }

  genes_keep <- intersect(colnames(datExpr_ref), colnames(datExpr_query))
  datExpr_ref <- datExpr_ref[, genes_keep]
  datExpr_query <- datExpr_query[, genes_keep]

  # set up multiExpr:
  setLabels <- c('ref', 'query')
  multiExpr <- list(
    ref = list(data = datExpr_ref),
    query = list(data = datExpr_query)
  )
  ref_modules <- GetModules(seurat_ref, wgcna_name = wgcna_name_ref)
  ref_modules <- ref_modules[genes_keep, ]
  ref_modules <- list(ref = ref_modules$module)



  # run the module preservation test:
  mp <- WGCNA::modulePreservation(
    multiExpr,
    ref_modules,
    referenceNetworks = 1,
    nPermutations = n_permutations,
    randomSeed = seed,
    quickCor = 0,
    parallelCalculation = parallel,
    ...
  )

  if (return_raw) {
    return(mp)
  }

  # get the stats obs and stats Z tables
  ref <- 1
  test <- 2
  statsObs <- cbind(
    mp$quality$observed[[ref]][[test]],
    mp$preservation$observed[[ref]][[test]]
  )
  statsZ <- cbind(
    mp$quality$Z[[ref]][[test]],
    mp$preservation$Z[[ref]][[test]]
  )

  # add stats to the seurat object:
  mod_pres <- list('obs' = statsObs, 'Z' = statsZ)
  seurat_obj <- SetModulePreservation(seurat_obj, mod_pres, name, wgcna_name)

  seurat_obj
}


#' ModulePreservationNetRep
#'
#' Computes module preservation statistics in a query dataset for a given reference dataset
#'
#' @param seurat_query A Seurat object
#' @param seurat_ref A Seurat object serving as the reference for the module preservation analysis
#' @param name The name to give the module preservation analysis.
#' @param n_permutations Number of permutations for the module preservation test.
#' @param n_threads Number of parallel threads to use for the module preservation test
#' @param TOM_use  The name of the hdWGCNA experiment containing the TOM that will be used for plotting
#' @param gene_mapping a dataframe containing gene name mappings between the query and the reference dataset. One column should have the gene name in the query dataset, and anotehr column should have the corresponding gene name in the reference dataset.
#' @param genome_ref_col the column name containing the gene names for the reference dataset
#' @param genome_query_col the column name containing the gene names for the query dataset
#' @param ref_power Soft power used to compute the reference adjacency matrix with WGCNA::adjacency. Default = 1
#' @param query_power Soft power used to compute the query adjacency matrix with WGCNA::adjacency. Default = 1
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_query@misc slot. Defaults to the active experiment.
#' @param wgcna_name_ref The name of the hdWGCNA experiment in the seurat_ref@misc slot. Defaults to the active experiment in seurat_ref.
#' @param ... Not currently used.
#' @details
#' ModulePreservationNetRep performs a statistical test to assess the preservation of co-expression modules identified
#' in one dataset in an independent dataset using the NetRep implementation rather than the WGCNA implementation.
#' This method is originally described by Langfelder et al. in the 2011 paper
#' "Is My Network Module Preserved and Reproducible?". The NetRep implementation of the module preservation test
#' is described by Ritchie et al. in the 2016 paper "A Scalable Permutation Approach Reveals Replication and Preservation Patterns of Networks Modules in Large Datasets"
#' This method can be used to assess biological differences in networks, as well as technical differences / reproducibility across different
#' batches. Note that the outputs of ModulePreservation and ModulePreservationNetRep are not identical. This function requires separate
#' installation of the NetRep package, which is not directly included as a dependency in hdWGCNA.
#'
#' @return The input seurat_query with the NetRep module preservation results stored under name in \code{seurat_query@misc[[wgcna_name]]}; retrieve with GetModulePreservation().
#' @export
ModulePreservationNetRep <- function(
  seurat_query,
  seurat_ref,
  name,
  n_permutations = 10000,
  n_threads = 8,
  gene_mapping = NULL,
  genome_ref_col = NULL,
  genome_query_col = NULL,
  TOM_use = NULL,
  ref_power = 1,
  query_power = 1,
  wgcna_name = NULL,
  wgcna_name_ref = NULL,
  ...
) {
  check_suggests('NetRep', feature = 'ModulePreservationNetRep')

  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_query@misc$active_wgcna
  }
  CheckWGCNAName(seurat_query, wgcna_name)
  if (is.null(TOM_use)) {
    TOM_use <- wgcna_name
  }

  if (is.null(wgcna_name_ref)) {
    wgcna_name_ref <- seurat_ref@misc$active_wgcna
  }
  CheckWGCNAName(seurat_ref, wgcna_name_ref)

  # get the modules from the reference dataset
  modules <- GetModules(seurat_ref, wgcna_name = wgcna_name_ref) %>% subset(module != 'grey')
  module_labels <- as.character(modules$module)
  genes_use_ref <- modules$gene_name
  names(module_labels) <- genes_use_ref

  # get datExpr for reference and query:
  datExpr_ref <- as.matrix(GetDatExpr(seurat_ref, wgcna_name_ref))

  # does this work? could be different genes.
  datExpr_query <- as.matrix(GetDatExpr(seurat_query, TOM_use))

  # get TOMs for reference and query:
  TOM_ref <- GetTOM(seurat_ref, wgcna_name_ref)
  TOM_query <- GetTOM(seurat_query, TOM_use)

  # change the gene names to match:
  if (!is.null(gene_mapping)) {
    gene_match <- match(colnames(datExpr_query), gene_mapping[, genome_query_col])
    gene_mapping <- na.omit(gene_mapping[gene_match, ])
    colnames(datExpr_query) <- gene_mapping[, genome_ref_col]
  }

  # genes to use in the query dataset
  genes_use_query <- genes_use_ref[genes_use_ref %in% colnames(TOM_query)]

  datExpr_ref <- datExpr_ref[, genes_use_ref]
  datExpr_query <- datExpr_query[, genes_use_query]
  TOM_ref <- TOM_ref[genes_use_ref, genes_use_ref]
  TOM_query <- TOM_query[genes_use_query, genes_use_query]

  # calculate adjacency:
  cor_mat_ref <- WGCNA::adjacency(datExpr_ref, type = 'unsigned', power = ref_power)
  cor_mat_query <- WGCNA::adjacency(datExpr_query, type = 'unsigned', power = query_power)

  # set up data lists:
  n_list <- list(TOM_ref, TOM_query)
  d_list <- list(datExpr_ref, datExpr_query)
  c_list <- list(cor_mat_ref, cor_mat_query)
  names(n_list) <- c('ref', 'query')
  names(d_list) <- c('ref', 'query')
  names(c_list) <- c('ref', 'query')

  # Assess the preservation of modules in the test dataset.
  preservation_test <- NetRep::modulePreservation(
    network = n_list,
    data = d_list,
    correlation = c_list,
    moduleAssignments = module_labels,
    discovery = 'ref',
    test = 'query',
    nPerm = n_permutations,
    nThreads = n_threads
  )

  # add the result to the seurat object
  seurat_query <- SetModulePreservation(
    seurat_query,
    preservation_test,
    name,
    wgcna_name
  )

  seurat_query
}


#' PlotModulePreservation
#'
#' Plotting function for Module Preservation statistics
#'
#' @param seurat_obj A Seurat object
#' @param name The name of the module preservation analysis to plot given in ModulePreservation
#' @param statistics Which module preservation statistics to plot? Choices are summary, rank, all, or a custom character vector of statistic names
#' @param plot_labels logical determining whether to plot the module labels
#' @param label_size the size of the module labels
#' @param mod_point_size the size of the points in each plot
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A ggplot object if a single statistic is plotted, otherwise a named list of ggplot objects (one per statistic).
#'
#' @details
#' This function creates a scatter plot showing the module preservation statistics for each module
#' compared to the size of the module (number of genes).
#'
#' @export
PlotModulePreservation <- function(
  seurat_obj,
  name,
  statistics = 'summary', # can be summary, rank, all, or a custom list
  plot_labels = TRUE,
  label_size = 4,
  mod_point_size = 4,
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get the module preservation stats:
  mod_pres <- GetModulePreservation(seurat_obj, name, wgcna_name)
  obs_df <- mod_pres$obs
  Z_df <- mod_pres$Z

  # get module colors:
  modules <- GetModules(seurat_obj, wgcna_name)
  module_colors <- modules %>%
    dplyr::select(c(module, color)) %>%
    dplyr::distinct()
  mods <- rownames(Z_df)
  mod_colors <- module_colors$color[match(mods, module_colors$module)]
  mod_colors <- ifelse(is.na(mod_colors), 'gold', mod_colors)

  # what are we going to plot?
  if (statistics == 'summary') {
    stat_list <- c('Zsummary.qual', 'Zsummary.pres')
  } else if (statistics == 'rank') {
    stat_list <- colnames(obs_df[, -1])[grepl('Rank', colnames(obs_df[, -1]))]
  } else if (statistics == 'all') {
    stat_list <- c(colnames(obs_df[, -1])[grepl('Rank', colnames(obs_df[, -1]))], colnames(Z_df[, -1]))
  } else {
    stat_list <- statistics
  }

  stat_list <- stat_list[stat_list != 'moduleSize']

  plot_list <- list()
  for (statistic in stat_list) {

    if (statistic %in% colnames(obs_df)) {
      values <- obs_df[, statistic]
    } else if (statistic %in% colnames(Z_df)) {
      values <- Z_df[, statistic]
    } else {
      stop('Invalid name for statistic.')
    }

    # setup plotting df
    plot_df <- data.frame(
      module = mods,
      color = mod_colors,
      value = values,
      size = Z_df$moduleSize
    )

    # don't include grey & gold:
    plot_df <- plot_df %>% subset(!(module %in% c('grey', 'gold')))

    if (grepl('Rank', statistic)) {
      cur_p <- plot_df %>% ggplot(aes(x = size, y = value, fill = module, color = module)) +
        geom_point(size = mod_point_size, pch = 21, color = 'black') +
        scale_y_reverse()
    } else {
      cur_p <- plot_df %>% ggplot(aes(x = size, y = value, fill = module, color = module)) +
        geom_rect(
          data = plot_df[1, ],
          aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 2), fill = 'grey75', alpha = 0.8, color = NA
        ) +
        geom_rect(
          data = plot_df[1, ],
          aes(xmin = 0, xmax = Inf, ymin = 2, ymax = 10), fill = 'grey92', alpha = 0.8, color = NA
        ) +
        geom_point(size = mod_point_size, pch = 21, color = 'black')
    }

    cur_p <- cur_p +
      scale_fill_manual(values = plot_df$color) +
      scale_color_manual(values = plot_df$color) +
      scale_x_continuous(trans = 'log10') +
      ylab(statistic) +
      xlab('Module Size') +
      ggtitle(statistic) +
      NoLegend() +
      theme(
        plot.title = element_text(hjust = 0.5)
      )


    if (plot_labels) {
      cur_p <- cur_p + ggrepel::geom_text_repel(label = plot_df$module, size = label_size, max.overlaps = Inf, color = 'black')
    }

    plot_list[[statistic]] <- cur_p
  }

  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  }

  plot_list
}


#' PlotModulePreservationLollipop
#'
#' Plots a ranked lollipop plot of co-expression modules based on the results of module preservation analysis.
#'
#' @return ggplot object containing the PlotModulePreservationLollipop
#'
#' @param seurat_obj A Seurat object
#' @param name The name to give the module preservation analysis.
#' @param features The name of the module preservation features to plot.
#' @param fdr logical indicating whether or not to plot FDR-corrected p-values
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' PlotModulePreservationLollipop generates a lollipop plot showing module preservation results. If the module
#' preservation test was performed using the WGCNA method, the statistic that will be shown is the Z-summary
#' preservation statistic. If the analysis was performed using NetRep, then the statistic that will be shown is
#' the FDR corrected averaged p-values from the module preservation permutation test.
#'
#' @export
PlotModulePreservationLollipop <- function(
  seurat_obj,
  name,
  features = NULL,
  fdr = TRUE,
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get the module info
  modules <- GetModules(seurat_obj, wgcna_name)
  mod_colors <- modules %>%
    dplyr::select(module, color) %>%
    dplyr::distinct()
  mod_cp <- mod_colors$color
  names(mod_cp) <- mod_colors$module

  # get the module preservation stats:
  mod_pres <- GetModulePreservation(seurat_obj, name, wgcna_name)

  # was this made with NetRep or WGCNA?
  if (all(names(mod_pres) %in% c('nulls', 'observed', 'p.values', 'nVarsPresent', 'propVarsPresent', 'totalSize', 'alternative'))) {
    netrep_used <- TRUE
  } else {
    netrep_used <- FALSE
  }

  if (netrep_used) {
    if (is.null(features)) {
      features <- 'average'
    }

    # check for valid features
    if (!features %in% c('average', colnames(mod_pres$p.value))) {
      stop(paste0('Invalid selection for features. Valid selections are: ', paste(c('average', colnames(mod_pres$p.value)), collapse = ', ')))
    }

    # pvals
    plot_df <- reshape2::melt(mod_pres$p.value)
    plot_df$Var2 <- paste0('pval.', plot_df$Var2)
    plot_df$type <- 'pval'
    plot_df1 <- plot_df

    # compute FDR:
    fdrs <- p.adjust(plot_df1$value, 'fdr')
    plot_df_fdrs <- plot_df1
    plot_df_fdrs$value <- fdrs
    plot_df_fdrs$type <- 'fdr'

    # observed stats
    plot_df <- reshape2::melt(mod_pres$observed)
    plot_df$type <- 'observation'
    plot_df2 <- plot_df

    plot_df <- rbind(plot_df1, plot_df_fdrs, plot_df2) %>%
      dplyr::rename(module = Var1, stat = Var2)

    # compute the average p-vals:
    if (features == 'average') {
      plot_title <- 'Summary'
      if (fdr) {
        plot_df <- plot_df %>%
          subset(type == 'fdr') %>%
          group_by(module) %>%
          summarise(value = mean(value))
        label <- bquote('-log'[10] ~ '(Avg. FDR)')
      } else {
        plot_df <- plot_df %>%
          subset(type == 'pval') %>%
          group_by(module) %>%
          summarise(value = mean(value))
        label <- bquote('-log'[10] ~ '(Avg. P-value)')
      }
    } else {
      plot_title <- features
      plot_df <- plot_df %>% subset(stat == paste0('pval.', features))
      if (fdr) {
        plot_df <- plot_df %>% subset(type == 'fdr')
        label <- bquote('-log'[10] ~ '(FDR)')
      } else {
        plot_df <- plot_df %>% subset(type == 'pval')
        label <- bquote('-log'[10] ~ '(P-value)')
      }
    }

    # add info about the module size (number of genes)
    mod_sizes <- mod_pres$nVarsPresent
    plot_df$mod_size <- mod_sizes[plot_df$module]

    # order by p-val
    plot_df <- plot_df %>% dplyr::arrange(desc(value))
    plot_df$module <- factor(as.character(plot_df$module), levels = as.character(plot_df$module))

    p <- plot_df %>%
      ggplot(aes(y = module, x = -log10(value), size = mod_size, color = module, fill = module)) +
      geom_rect(
        data = plot_df[1, ],
        aes(xmin = -Inf, ymax = Inf, ymin = -Inf, xmax = -log10(0.05)), fill = 'grey90', alpha = 0.8, color = NA
      ) +
      # geom_vline(xintercept=2, linewidth=0.5, linetype='dashed') +
      geom_segment(aes(y = module, yend = module, x = 0, xend = -log10(value)), linewidth = 0.5, alpha = 0.5) +
      geom_point(shape = 21, color = 'black') +
      scale_color_manual(values = mod_cp, guide = 'none') +
      scale_fill_manual(values = mod_cp, guide = 'none') +
      ylab('') +
      labs(size = bquote('N'[genes])) +
      xlab(label) +
      theme(
        panel.border = element_rect(linewidth = 1, color = 'black', fill = NA),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = 'bold')
      ) +
      ggtitle(plot_title)
  } else {
    mod_pres <- mod_pres$Z

    if (is.null(features)) {
      features <- 'Zsummary.pres'
    }

    # check for valid features
    if (!features %in% colnames(mod_pres)) {
      stop(paste0('Invalid selection for features. Valid selections are: ', paste(c('average', colnames(mod_pres), collapse = ', '))))
    }

    plot_df <- mod_pres[, c('moduleSize', features)]
    colnames(plot_df)[2] <- 'value'
    plot_df$module <- rownames(plot_df)
    plot_df <- plot_df %>% subset(!module %in% c('gold', 'grey'))
    plot_df <- plot_df %>% dplyr::arrange(value, descending = TRUE)
    plot_df$module <- factor(as.character(plot_df$module), levels = as.character(plot_df$module))

    p <- plot_df %>%
      ggplot(aes(y = module, x = value, size = moduleSize, color = module, fill = module)) +
      geom_rect(
        data = plot_df[1, ],
        aes(xmin = -Inf, ymax = Inf, ymin = -Inf, xmax = 2), fill = 'grey75', alpha = 0.8, color = NA
      ) +
      geom_rect(
        data = plot_df[1, ],
        aes(ymin = -Inf, ymax = Inf, xmin = 2, xmax = 10), fill = 'grey92', alpha = 0.8, color = NA
      ) +
      geom_segment(aes(y = module, yend = module, x = 0, xend = value), linewidth = 0.5, alpha = 0.5) +
      geom_point(shape = 21, color = 'black') +
      scale_color_manual(values = mod_cp, guide = 'none') +
      scale_fill_manual(values = mod_cp, guide = 'none') +
      ylab('') +
      xlab('') +
      ggtitle(features) +
      theme(
        panel.border = element_rect(linewidth = 1, color = 'black', fill = NA),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = 'bold')
      ) +
      RotatedAxis()
  }

  # return the plot
  p
}
