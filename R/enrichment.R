#' RunEnrichr
#'
#' Run Enrichr gene set enrichment tests on hdWGCNA modules
#'
#' @param seurat_obj A Seurat object
#' @param dbs character vector of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wait logical indicating whether or not to wait some time between sending requests to the EnrichR server.
#' @param wait_time the number of seconds to wait between sending requests to the EnrichR server. Value must be less than 60.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Currently unused.
#' @return The Seurat object with the Enrichr results stored in \code{seurat_obj@misc[[wgcna_name]]}; retrieve them with \code{GetEnrichrTable}.
#' @keywords scRNA-seq
#' @export
#' RunEnrichr
RunEnrichr <- function(
  seurat_obj,
  dbs = c('GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021'),
  max_genes = 100,
  wait = TRUE,
  wait_time = 5,
  wgcna_name = NULL, ...
) {
  check_suggests('enrichR', feature = 'RunEnrichr')
  old_options <- enrichr_options()
  on.exit(options(old_options), add = TRUE)

  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # check wait_time
  if (!is.numeric(wait_time)) {
    stop(paste0('wait_time must be a numeric.'))
  }
  if (wait_time > 60 || wait_time < 1) {
    stop(paste0('Invalid value selected for wait_time, must be greater than 0 and less than 60.'))
  }

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)

  # exclude grey
  mods <- mods[mods != 'grey']

  # hub genes
  hub_df <- GetHubGenes(seurat_obj, n_hubs = max_genes, wgcna_name = wgcna_name)

  # run EnrichR for loop:
  combined_output <- data.frame()
  for (i in seq_along(mods)) {
    cur_mod <- mods[i]
    cur_info <- subset(hub_df, module == cur_mod)
    cur_genes <- as.character(cur_info$gene_name)

    # run the enrichment test
    enriched <- enrichR::enrichr(cur_genes, dbs)

    if (wait) {
      Sys.sleep(wait_time)
    }

    # collapse into one db
    for (db in names(enriched)) {
      cur_df <- enriched[[db]]
      if (nrow(cur_df) > 1) {
        cur_df$db <- db
        cur_df$module <- cur_mod
        combined_output <- rbind(combined_output, cur_df)
      }
    }
  }

  # add GO term table to seurat object
  seurat_obj <- SetEnrichrTable(seurat_obj, combined_output, wgcna_name)
  seurat_obj
}


#' EnrichrBarPlot
#'
#' Generates bar plots from Enrichr data to visualize enriched terms for hdWGCNA modules in a Seurat object.
#' The function outputs a PDF file for each module, with separate bar plots for each database.
#'
#' @param seurat_obj A Seurat object
#' @param outdir A string specifying the directory where the output PDF files will be saved. Default is "enrichr_plots".
#' @param n_terms An integer indicating the number of top enriched terms to include in each bar plot. Default is 25.
#' @param p_cutoff A numeric value specifying the significance threshold for including terms (p-value or adjusted p-value).
#'                 Only terms with p-values below this threshold are plotted. Default is 0.05.
#' @param p_adj Logical indicating whether to use the adjusted p-value (default: TRUE) or raw p-value for filtering terms.
#' @param plot_size A numeric vector of length 2 specifying the width and height of the output PDF files in inches. Default is c(6, 15).
#' @param logscale Logical specifying whether to log-transform the enrichment scores before plotting. Default is FALSE.
#' @param plot_bar_color A string specifying the color of the bars in the bar plots. If NULL (default), bars are colored according
#'                       to the module's assigned color.
#' @param plot_text_color A string specifying the color of the text labels on the bar plots. If NULL (default), the color is
#'                        automatically determined based on the bar color.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#'
#' @details
#' This function processes the Enrichr output stored in a Seurat object, filters enriched terms by significance, and generates
#' bar plots for each WGCNA module. Separate plots are created for each database included in the Enrichr results. The top
#' enriched terms for each module and database are ordered by their combined enrichment score. If there are ties in the enrichment
#' scores, the function uses all tied terms.
#'
#' The bar plot text and bar colors can be customized, and plots can be saved with log-transformed enrichment scores if specified.
#' Text wrapping is applied to long term names for better readability.
#'
#' @return
#' Generates PDF files in the specified output directory, with one file per WGCNA module. Each file contains bar plots for
#' all databases associated with that module.
#'
#' @export
EnrichrBarPlot <- function(
  seurat_obj,
  outdir = 'enrichr_plots',
  n_terms = 25,
  p_cutoff = 0.05,
  p_adj = TRUE,
  plot_size = c(6, 15),
  logscale = FALSE,
  plot_bar_color = NULL,
  plot_text_color = NULL,
  wgcna_name = NULL
) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get Enrichr table
  enrichr_df <- GetEnrichrTable(seurat_obj, wgcna_name)
  dbs <- as.character(unique(enrichr_df$db))

  # subset based on significance level:
  if (p_adj) {
    enrichr_df <- subset(enrichr_df, Adjusted.P.value <= p_cutoff)
  } else {
    enrichr_df <- subset(enrichr_df, P.value <= p_cutoff)
  }

  # helper function to wrap text
  wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = '\n'), USE.NAMES = FALSE)
  }

  # make output dir if it doesn't exist:
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }

  # loop through modules:
  for (i in seq_along(mods)) {
    cur_mod <- mods[i]
    cur_terms <- subset(enrichr_df, module == cur_mod)
    print(cur_mod)

    # get color for this module:
    cur_color <- modules %>%
      subset(module == cur_mod) %>%
      .$color %>%
      unique() %>%
      as.character()
    if (!is.null(plot_bar_color)) {
      cur_color <- plot_bar_color
    }

    # skip if there are not any terms for this module:
    if (nrow(cur_terms) == 0) {
      next
    }
    cur_terms$wrap <- wrapText(cur_terms$Term, 45)

    # plot top n_terms as barplot
    plot_list <- list()
    for (cur_db in dbs) {
      plot_df <- subset(cur_terms, db == cur_db) %>%
        slice_max(order_by = Combined.Score, n = n_terms)

      # text color:
      if (is.null(plot_text_color)) {
        if (cur_color == 'black') {
          text_color <- 'grey'
        } else {
          text_color <- 'black'
        }
      } else {
        text_color <- plot_text_color
      }


      # logscale?
      if (logscale) {
        plot_df$Combined.Score <- log(plot_df$Combined.Score)
        lab <- 'Enrichment log(combined score)'
        x <- 0.2
      } else {
        lab <- 'Enrichment (combined score)'
        x <- 5
      }

      # make bar plot:
      plot_list[[cur_db]] <- ggplot(plot_df, aes(x = Combined.Score, y = reorder(wrap, Combined.Score))) +
        geom_bar(stat = 'identity', position = 'identity', color = 'white', fill = cur_color) +
        geom_text(aes(label = wrap), x = x, color = text_color, size = 3.5, hjust = 'left') +
        scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
        xlab(lab) +
        ylab('') +
        ggtitle(cur_db) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.line.y = element_blank()
        )
    }

    # make pdfs in output dir
    pdf(paste0(outdir, '/', cur_mod, '.pdf'), width = plot_size[1], height = plot_size[2])
    for (plot in plot_list) {
      print(plot)
    }
    dev.off()
  }
}


#' EnrichrDotPlot
#'
#' Generate a dot plot visualizing enrichment results from Enrichr for hdWGCNA modules.
#'
#' This function creates dot plots from Enrichr results associated with hdWGCNA modules.
#' Each module is represented by its most enriched terms from the specified Enrichr database.
#' The size of the dots indicates the enrichment score, and the color indicates the statistical significance (-log10 transformed p-value).
#'
#' @param seurat_obj A Seurat object
#' @param database A character string specifying the name of the Enrichr database to use (e.g., "GO_Biological_Process_2021").
#' @param mods A character vector specifying the names of modules to include in the plot.
#'        If `mods = "all"` (default), all modules except the "grey" module are included.
#' @param n_terms An integer specifying the number of top enriched terms to plot for each module (default = 3).
#' @param p_cutoff A numeric value specifying the p-value threshold for filtering enriched terms (default = 0.05).
#' @param p_adj A logical value indicating whether to use adjusted p-values (`TRUE`, default) or raw p-values (`FALSE`) for filtering.
#' @param break_ties A logical value indicating whether to randomly select among tied terms to enforce `n_terms` (default = `TRUE`).
#' @param term_size A numeric value specifying the font size of the enriched terms displayed on the y-axis (default = 10).
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' - The function first retrieves WGCNA module and Enrichr data from the specified Seurat object.
#' - Modules are filtered based on the `mods` parameter, and enriched terms are filtered by significance using `p_cutoff` and `p_adj`.
#' - The top `n_terms` terms for each module are selected based on the Combined Score. If ties occur, the `break_ties` parameter determines how they are resolved.
#' - A dot plot is generated where each dot represents an enriched term, its size corresponds to the Combined Score (log-transformed),
#'   and its color indicates the significance (-log10 transformed p-value).
#'
#' @return A ggplot2 object representing the dot plot of enriched terms for the specified modules and database.
#' @export
#'
EnrichrDotPlot <- function(
  seurat_obj,
  database,
  mods = 'all',
  n_terms = 3,
  p_cutoff = 0.05,
  p_adj = TRUE,
  break_ties = TRUE,
  term_size = 10,
  wgcna_name = NULL
) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)

  # using all modules?
  if (mods == 'all') {
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']
  }

  # get Enrichr table
  enrichr_df <- GetEnrichrTable(seurat_obj, wgcna_name)

  # --- ✅ CHECK 1: validate that the requested database exists ---
  if (!'db' %in% colnames(enrichr_df)) {
    stop("The Enrichr table does not contain a 'db' column. Please check the enrichment results.")
  }

  available_dbs <- unique(enrichr_df$db)
  if (!(database %in% available_dbs)) {
    stop(
      paste0(
        "Database '", database, "' not found in enrichment results.\n",
        'Available databases: ', paste(available_dbs, collapse = ', ')
      )
    )
  }

  # subset based on significance level
  if (p_adj) {
    enrichr_df <- subset(enrichr_df, Adjusted.P.value <= p_cutoff)
  } else {
    enrichr_df <- subset(enrichr_df, P.value <= p_cutoff)
  }

  # add color to enrich_table
  mod_colors <- dplyr::select(modules, c(module, color)) %>% dplyr::distinct()
  enrichr_df$color <- mod_colors[match(enrichr_df$module, mod_colors$module), 'color']

  # helper function to wrap text
  wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = '\n'), USE.NAMES = FALSE)
  }

  # get top terms
  top_terms <- enrichr_df %>%
    subset(db == database & module %in% mods) %>%
    dplyr::group_by(module) %>%
    dplyr::slice_max(order_by = Combined.Score, n = n_terms) %>%
    .$Term

  plot_df <- subset(enrichr_df, Term %in% top_terms & db == database & module %in% mods)

  # --- ✅ CHECK 2: warn if no terms or missing modules ---
  if (nrow(plot_df) == 0) {
    warning(
      paste0(
        "No enriched terms to plot for database '", database,
        "' with p_cutoff = ", p_cutoff,
        ' and selected modules: ', paste(mods, collapse = ', ')
      )
    )
    return(NULL)
  }

  # identify modules missing any significant terms
  plotted_modules <- unique(plot_df$module)
  missing_modules <- setdiff(mods, plotted_modules)
  if (length(missing_modules) > 0) {
    warning(
      paste0(
        "No enriched terms found for the following modules in database '", database, "': ",
        paste(missing_modules, collapse = ', ')
      )
    )
  }

  # randomly break ties if requested
  if (break_ties) {
    plot_df <- do.call(
      rbind,
      lapply(
        plot_df %>% dplyr::group_by(module) %>% dplyr::group_split(),
        function(x) {
          if (nrow(x) > n_terms) x <- x[sample(seq_len(nrow(x)), n_terms), ]
          return(x)
        }
      )
    )
  }

  plot_df <- plot_df %>%
    dplyr::mutate(Term = stringr::str_replace(Term, ' \\s*\\([^\\)]+\\)', ''))
  plot_df$Term <- wrapText(plot_df$Term, 45)

  # set module factor order
  plot_df$module <- factor(
    as.character(plot_df$module),
    levels = levels(modules$module)
  )
  plot_df <- dplyr::arrange(plot_df, module)

  # set Terms factor
  plot_df$Term <- factor(
    as.character(plot_df$Term),
    levels = unique(as.character(plot_df$Term))
  )

  # p-value column
  plot_df$p <- if (p_adj) plot_df$Adjusted.P.value else plot_df$P.value
  plot_df$p <- ifelse(is.na(plot_df$p), 1, plot_df$p)

  max_p <- quantile(-log(plot_df$p), 0.95)
  plot_df$logp <- pmin(-log(plot_df$p), max_p)

  p <- plot_df %>%
    ggplot(aes(x = module, y = Term, size = log10(Combined.Score), color = logp)) +
    geom_point() +
    Seurat::RotatedAxis() +
    ylab('') +
    xlab('') +
    labs(
      color = bquote('-log'[10] ~ '(P)'),
      size  = bquote('log'[10] ~ '(Enrich)')
    ) +
    scale_y_discrete(limits = rev) +
    ggtitle(database) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_text(size = term_size),
      panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1),
      panel.grid = element_line(linewidth = 0.25, color = 'lightgrey')
    )

  return(p)
}


#' OverlapModulesDEGs
#'
#' Performs Fisher's Exact Test for overlap between DEGs and hdWGCNA modules.
#'
#' @param seurat_obj A Seurat object
#' @param deg_df DEG table formatted like the output from Seurat's FindMarkers
#' @param fc_cutoff log fold change cutoff for DEGs to be included in the overlap test
#' @param group_col the name of the column in deg_df containing the cell grouping information
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Currently unused.
#' @return A data.frame with one row per module and cell group, containing the module color, odds ratio, p-value, FDR, significance stars, Jaccard index, and overlap size.
#' @keywords scRNA-seq
#' @export
#' OverlapModulesDEGs
OverlapModulesDEGs <- function(
  seurat_obj,
  deg_df,
  fc_cutoff = 0.5,
  group_col = 'cluster',
  wgcna_name = NULL,
  ...
) {
  check_suggests('GeneOverlap', feature = 'OverlapModulesDEGs')

  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  deg_df$group <- deg_df[, group_col]
  cell_groups <- deg_df$group %>% unique()

  # get modules,
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # subset deg_df by fold change cutoff:
  if (fc_cutoff >= 0) {
    deg_df <- subset(deg_df, avg_log2FC >= fc_cutoff)
  } else {
    deg_df <- subset(deg_df, avg_log2FC <= fc_cutoff)

    # reverse the sign of the remaining fold changes
    deg_df$avg_log2FC <- -1 * deg_df$avg_log2FC
    fc_cutoff <- -1 * fc_cutoff
  }

  # size of genome based on # genes in Seurat object:
  genome.size <- nrow(seurat_obj)

  # run overlaps between module gene lists and DEG lists:
  overlap_df <- do.call(rbind, lapply(mods, function(cur_mod) {
    cur_module_genes <- modules %>%
      subset(module == cur_mod) %>%
      .$gene_name
    cur_overlap_df <- do.call(rbind, lapply(cell_groups, function(cur_group) {
      cur_DEGs <- deg_df %>%
        subset(group == cur_group & p_val_adj <= 0.05 & avg_log2FC > fc_cutoff) %>%
        .$gene
      cur_overlap <- GeneOverlap::testGeneOverlap(GeneOverlap::newGeneOverlap(
        cur_module_genes,
        cur_DEGs,
        genome.size = genome.size
      ))
      c(cur_overlap@odds.ratio, cur_overlap@pval, cur_overlap@Jaccard, length(cur_overlap@intersection))
    })) %>% as.data.frame()
    colnames(cur_overlap_df) <- c('odds_ratio', 'pval', 'Jaccard', 'size_intersection')
    cur_overlap_df$module <- cur_mod
    cur_overlap_df$group <- cell_groups

    # module color:
    cur_overlap_df$color <- modules %>%
      subset(module == cur_mod) %>%
      .$color %>%
      unique()
    cur_overlap_df
  }))


  # adjust for multiple comparisons:
  overlap_df$fdr <- p.adjust(overlap_df$pval, method = 'fdr')

  # significance level:
  overlap_df$Significance <- gtools::stars.pval(overlap_df$fdr)
  overlap_df$Significance <- ifelse(
    overlap_df$Significance == '.', '',
    overlap_df$Significance
  )

  # set factor levels for modules:
  overlap_df$module <- factor(overlap_df$module, levels = mods)

  # re-arrange columns:
  overlap_df <- overlap_df %>% dplyr::select(c(module, group, color, odds_ratio, pval, fdr, Significance, Jaccard, size_intersection))

  overlap_df
}


#' OverlapDotPlot
#'
#' Makes barplots from Enrichr data
#'
#' @param overlap_df the Module/DEG overlap table from OverlapModulesDEGs
#' @param plot_var the name of the overlap statistic to plot
#' @param logscale logical controlling whether to plot the result on a log scale, useful for odds ratio
#' @param neglog logical controlling wehether to plot the result as a negative log, useful for p-value / FDR
#' @param plot_significance logical controlling whether to plot the significance levels on top of the dots
#' @param ... Currently unused.
#' @return A ggplot object.
#' @keywords scRNA-seq
#' @export
OverlapDotPlot <- function(
  overlap_df, plot_var = 'odds_ratio',
  logscale = TRUE,
  neglog = FALSE,
  plot_significance = TRUE,
  ...
) {
  label <- plot_var
  if (logscale) {
    overlap_df[[plot_var]] <- log(overlap_df[[plot_var]])
    label <- paste0('log(', plot_var, ')')
  }
  if (neglog) {
    overlap_df[[plot_var]] <- -1 * overlap_df[[plot_var]]
    label <- paste0('-', label)
  }

  p <- overlap_df %>% ggplot(aes(x = module, y = group)) +
    geom_point(
      aes(
        size = get(plot_var),
        alpha = get(plot_var)
      ),
      color = overlap_df$color
    ) +
    RotatedAxis() +
    ylab('') +
    xlab('') +
    labs(size = label, alpha = label) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1)
    )

  # plot significance level?
  if (plot_significance) {
    p <- p + geom_text(aes(label = Significance))
  }

  p
}


#' OverlapBarPlot
#'
#' Plots the results from OverlapModulesDEGs as a bar plot
#'
#' @param overlap_df the Module/DEG overlap table from OverlapModulesDEGs
#' @param plot_var the name of the overlap statistic to plot
#' @param logscale logical controlling whether to plot the result on a log scale, useful for odds ratio
#' @param neglog logical controlling wehether to plot the result as a negative log, useful for p-value / FDR
#' @param label_size the size of the module labels in the bar plot
#' @param ... Currently unused.
#' @return A named list of ggplot objects, one bar plot per cell group.
#' @keywords scRNA-seq
#' @export
OverlapBarPlot <- function(
  overlap_df,
  plot_var = 'odds_ratio',
  logscale = FALSE, neglog = FALSE,
  label_size = 2,
  ...
) {
  label <- plot_var

  if (plot_var == 'odds_ratio') {
    yint <- 1
  } else if (plot_var == 'fdr') {
    yint <- 0.05
  }

  if (logscale) {
    overlap_df[[plot_var]] <- log(overlap_df[[plot_var]])
    label <- paste0('log(', plot_var, ')')
    yint <- log(yint)
  }
  if (neglog) {
    overlap_df[[plot_var]] <- -1 * log(overlap_df[[plot_var]])
    label <- paste0('-log(', label, ')')
    yint <- -1 * log(yint)
  }

  groups <- overlap_df$group %>%
    as.character() %>%
    unique()

  plot_list <- list()
  for (cur_group in groups) {
    cur_df <- overlap_df %>%
      subset(group == cur_group)

    p <- cur_df %>%
      ggplot(aes(x = reorder(module, get(plot_var)), y = get(plot_var))) +
      geom_bar(stat = 'identity', fill = cur_df$color) +
      coord_flip() +
      xlab('') +
      ylab(label) +
      ggtitle(cur_group) +
      theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )

    if (plot_var == 'fdr' || plot_var == 'odds_ratio') {
      p <- p + geom_hline(yintercept = yint, linetype = 'dashed', color = 'gray')
    }

    # add the labels:
    p <- p +
      geom_text(
        aes(label = module, x = module, y = get(plot_var)),
        color = 'black', size = label_size, hjust = 'inward'
      )

    plot_list[[cur_group]] <- p
  }

  plot_list
}

# enrichR sets its server address in .onAttach(), so calling it through the
# namespace without library(enrichR) fails with 'Must specify at least one of
# url or handle'. Fill in the defaults when they are missing; the caller
# restores the previous options on exit.
enrichr_options <- function() {
  if (!is.null(getOption('enrichR.base.address'))) {
    return(list())
  }
  base <- 'https://maayanlab.cloud/'
  options(
    enrichR.sites.base.address = base,
    enrichR.base.address = paste0(base, 'Enrichr/'),
    speedrichr.base.address = paste0(base, 'speedrichr/api/'),
    enrichR.live = TRUE,
    enrichR.quiet = TRUE,
    modEnrichR.use = TRUE,
    enrichR.sites = c('Enrichr', 'FlyEnrichr', 'WormEnrichr', 'YeastEnrichr', 'FishEnrichr', 'OxEnrichr')
  )
}
