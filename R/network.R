#' TestSoftPowers
#'
#' Compute the scale-free topology model fit for different soft power thresholds
#'
#' @param seurat_obj A Seurat object
#' @param powers numeric vector specifying soft powers to test
#' @param networkType The type of network to use for network analysis. Options are "signed" (default), "unsigned", or "signed hybrid". This should be consistent with the network chosen for ConstructNetwork
#' @param corFnc Correlation function for the gene-gene correlation adjacency matrix.
#' @param wgcna_name name of the WGCNA experiment
#' @param ... additional parameters passed to WGCNA::pickSoftThreshold
#'
#' @details
#' TestSoftPowers aims to find the optimal value of the "soft-power threshold" for co-expression network analysis.
#' For varying values of this soft_power parameter, this function assesses the scale-free topology model fit of the resulting network
#' and the function also calculates other network statistics such as the connectivity. The user may also vary the type of network used and the correlation function.
#' For downstream functions, the lowest soft_power that reaches a fit of 0.8 is selected for network construction, but the user may select their own values.
#' @export
TestSoftPowers <- function(
  seurat_obj,
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  networkType = 'signed',
  corFnc = 'bicor',
  wgcna_name = NULL,
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get datExpr
  datExpr <- GetDatExpr(seurat_obj, wgcna_name)

  # Call the network topology analysis function
  powerTable <- list(
    data = WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      verbose = 100,
      networkType = networkType,
      corFnc = corFnc,
      ...
    )[[2]]
  )

  # set the power table in Seurat object:
  seurat_obj <- SetPowerTable(seurat_obj, powerTable$data, wgcna_name = wgcna_name)
  seurat_obj
}


#' TestSoftPowersConsensus
#'
#' Compute the scale-free topology model fit for different soft power thresholds separately for each input dataset
#'
#' @param seurat_obj A Seurat object
#' @param powers numeric vector specifying soft powers to test
#' @param use_metacells logical flag for whether to use the metacell expression matrix
#' @param networkType The type of network to use for network analysis. Options are "signed" (default), "unsigned", or "signed hybrid". This should be consistent with the network chosen for ConstructNetwork
#' @param corFnc Correlation function for the gene-gene correlation adjacency matrix.
#' @param setDatExpr logical flag indicating whether to run setDatExpr.
#' @param group.by A string containing the name of a column in the Seurat object with cell groups (clusters, cell types, etc). If NULL (default), hdWGCNA uses the Seurat Idents as the group.
#' @param group_name A string containing a group present in the provided group.by column or in the Seurat Idents. A character vector can be provided to select multiple groups at a time.
#' @param multi.group.by A string containing the name of a column in the Seurat object with groups for consensus WGCNA (dataset, sample, condition, etc)
#' @param multi_groups A character vecrtor containing the names of groups to select
#' @param ... additional parameters passed to SetDatExpr
#' @keywords scRNA-seq
#' @export
#' @examples
#' # TestSoftPowers(pbmc)
TestSoftPowersConsensus <- function(
  seurat_obj,
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  use_metacells = TRUE,
  networkType = 'signed',
  corFnc = 'bicor',
  setDatExpr = FALSE,
  group.by = NULL, group_name = NULL,
  multi.group.by = NULL,
  multi_groups = NULL,
  wgcna_name = NULL,
  ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # add multiExpr if not already added:
  if (!('multiExpr' %in% names(seurat_obj@misc[[wgcna_name]])) || setDatExpr == TRUE) {
    # set
    seurat_obj <- SetMultiExpr(
      seurat_obj,
      group_name = group_name,
      group.by = group.by,
      multi.group.by = multi.group.by,
      multi_groups = multi_groups,
      wgcna_name = wgcna_name
    )
  }

  multiExpr <- GetMultiExpr(seurat_obj, wgcna_name)

  # pick soft thresh for each consensus group:
  powerTables <- list()
  for (i in seq_along(multiExpr)) {
    cur_group <- names(multiExpr)[i]
    print(cur_group)

    # Call the network topology analysis function for each set in turn
    powerTable <- list(
      data = WGCNA::pickSoftThreshold(
        multiExpr[[cur_group]]$data,
        powerVector = powers,
        verbose = 100,
        networkType = networkType,
        corFnc = corFnc,
        ...
      )[[2]]
    )
    powerTable$data$group <- cur_group
    powerTables[[cur_group]] <- powerTable$data
  }

  # merge the power tables
  powerTable <- do.call(rbind, powerTables)

  # set the power table in Seurat object:
  seurat_obj <- SetPowerTable(seurat_obj, powerTable, wgcna_name = wgcna_name)
  seurat_obj
}


#' PlotSoftPowers
#'
#' Plot Soft Power Threshold results
#'
#' @param seurat_obj A Seurat object
#' @param selected_power power to highlight in the plots
#' @param point_size the size of the points in the plot
#' @param text_size the size of the text in the plot
#' @param plot_connectivity logical indicating whether to plot the connectivity in addition to the scale free topplogy fit.
#' @param wgcna_name The name of the WGCNA experiment in seurat_obj
#' @keywords scRNA-seq
#' @export
PlotSoftPowers <- function(
  seurat_obj,
  selected_power = NULL,
  point_size = 5,
  text_size = 3,
  plot_connectivity = TRUE,
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  pt <- GetPowerTable(seurat_obj, wgcna_name)

  if ('group' %in% colnames(pt)) {
    # select soft power for each group:
    power_tables <- pt %>% dplyr::group_split(group)
    soft_powers <- sapply(power_tables, function(power_table) {
      power_table %>%
        subset(SFT.R.sq >= 0.8) %>%
        .$Power %>%
        min()
    })
  } else {
    # get soft power:
    if (is.null(selected_power)) {
      soft_power <- pt %>%
        subset(SFT.R.sq >= 0.8) %>%
        .$Power %>%
        min()
    } else {
      soft_power <- selected_power
    }
    soft_powers <- NULL
    power_tables <- list('power' = pt)
  }

  plot_list <- list()

  for (i in seq_along(power_tables)) {
    pt <- power_tables[[i]]

    if (!is.null(soft_powers)) {
      soft_power <- soft_powers[i]
      print(i)
      print(soft_power)

      pt <- pt %>% dplyr::select(-group)
    }

    print(head(pt))

    # get other params:
    sft_r <- as.numeric(pt[pt$Power == soft_power, 'SFT.R.sq'])
    mean_k <- as.numeric(pt[pt$Power == soft_power, 'mean.k.'])
    median_k <- as.numeric(pt[pt$Power == soft_power, 'median.k.'])
    max_k <- as.numeric(pt[pt$Power == soft_power, 'max.k.'])

    # set color of text
    pt$text_color <- ifelse(
      pt$Power == soft_power, 'white', 'black'
    )

    # plot for soft power thresh:
    p1 <- pt %>% ggplot(aes(x = Power, y = SFT.R.sq)) +
      geom_rect(
        data = pt[1, ],
        aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.8), fill = 'grey80', alpha = 0.8, color = NA
      ) +
      geom_hline(yintercept = sft_r, linetype = 'dashed') +
      geom_vline(xintercept = soft_power, linetype = 'dashed') +
      geom_point(
        data = pt[pt$Power == soft_power, c('Power', 'SFT.R.sq')],
        aes(x = Power, y = SFT.R.sq),
        inherit.aes = FALSE,
        color = 'black',
        size = point_size
      ) +
      geom_text(label = pt$Power, color = pt$text_color, size = text_size) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      ylab('Scale-free Topology Model Fit') +
      xlab('Soft Power Threshold') +
      theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1)
      )

    if (plot_connectivity) {
      # plot for mean connectivity:
      p2 <- pt %>% ggplot(aes(x = Power, y = mean.k.)) +
        geom_hline(yintercept = mean_k, linetype = 'dashed') +
        geom_vline(xintercept = soft_power, linetype = 'dashed') +
        geom_point(
          data = pt[pt$Power == soft_power, c('Power', 'mean.k.')],
          aes(x = Power, y = mean.k.),
          inherit.aes = FALSE,
          color = 'black',
          size = point_size
        ) +
        geom_text(label = pt$Power, color = pt$text_color, size = text_size) +
        scale_y_continuous(labels = scales::comma) +
        ylab('Mean Connectivity') +
        xlab('Soft Power Threshold') +
        theme(
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1)
        )

      # plot for medianan connectivity:
      p3 <- pt %>% ggplot(aes(x = Power, y = median.k.)) +
        geom_hline(yintercept = median_k, linetype = 'dashed') +
        geom_vline(xintercept = soft_power, linetype = 'dashed') +
        geom_point(
          data = pt[pt$Power == soft_power, c('Power', 'median.k.')],
          aes(x = Power, y = median.k.),
          inherit.aes = FALSE,
          color = 'black',
          size = point_size
        ) +
        geom_text(label = pt$Power, color = pt$text_color, size = text_size) +
        scale_y_continuous(labels = scales::comma) +
        ylab('Median Connectivity') +
        xlab('Soft Power Threshold') +
        theme(
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1)
        )

      # plot for mean connectivity:
      p4 <- pt %>% ggplot(aes(x = Power, y = max.k.)) +
        geom_hline(yintercept = max_k, linetype = 'dashed') +
        geom_vline(xintercept = soft_power, linetype = 'dashed') +
        geom_point(
          data = pt[pt$Power == soft_power, c('Power', 'max.k.')],
          aes(x = Power, y = max.k.),
          inherit.aes = FALSE,
          color = 'black',
          size = point_size
        ) +
        geom_text(label = pt$Power, color = pt$text_color, size = text_size) +
        scale_y_continuous(labels = scales::comma) +
        ylab('Max Connectivity') +
        xlab('Soft Power Threshold') +
        theme(
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1)
        )


      plot_list[[i]] <- list(p1, p2, p3, p4)
    } else {
      plot_list[[i]] <- p1
    }
  }

  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  }
  plot_list
}


#' ConstructNetwork
#'
#' @description Constructs a co-expression network and groups genes into modules
#' given a Seurat object that has been prepared for network analysis.
#'
#' @return seurat_obj with the co-expression network and gene modules computed for the selected wgcna experiment
#'
#' @param seurat_obj A Seurat object
#' @param soft_power the soft power used for network construction. Automatically selected by default.
#' @param min_power the smallest soft power to be selected if soft_power=NULL
#' @param tom_outdir path to the directory where the TOM will be written
#' @param tom_name prefix name given to the TOM output file
#' @param consensus flag indicating whether or not to perform Consensus network analysis
#' @param wgcna_name name of the WGCNA experiment
#' @param ... additional parameters passed to blockwiseConsensusModules
#'
#' @details
#' ConstructNetwork builds a co-expression network and identifies clusters of
#' highly co-expressed genes (modules) from the metacell or metaspot
#' expression matrix stored in the Seurat object. Before running this function,
#' the following functions must be run on the input Seurat object:
#'
#' 1. SetupForWGCNA
#' 2. MetacellsByGroups or MetaspotsByGroups, and NormalizeMetacells
#' 3. SetDatExpr or SetMultiExpr
#' 4. TestSoftPowers or TestSoftPowersConsensus
#'
#' This function can also be used to perform consensus network analysis if consensus=TRUE
#' is selected. ConstructNetwork calls the WGCNA function blockwiseConsensusModules
#' to compute the adjacency matrix, topological overlap matrix, and to run the
#' Dynamic Tree Cut algorithm to identify gene modules. blockwiseConsensusModules
#' has numerous parameters but here we have selected default parameters that we
#' have found to provide reasonable results on a variety of single-cell and
#' spatial transcriptomic datasets.
#'
#' @export
ConstructNetwork <- function(
  seurat_obj, soft_power = NULL, min_power = 3,
  tom_outdir = 'TOM',
  tom_name = NULL,
  consensus = FALSE,
  overwrite_tom = FALSE,
  store_tom_in_seurat = FALSE,
  wgcna_name = NULL,
  blocks = NULL, maxBlockSize = 30000, randomSeed = 12345, corType = 'pearson',
  consensusQuantile = 0.3, networkType = 'signed', TOMType = 'signed',
  TOMDenom = 'min', scaleTOMs = TRUE, calibrationQuantile = 0.95,
  sampleForCalibration = TRUE, sampleForCalibrationFactor = 1000,
  useDiskCache = TRUE, chunkSize = NULL,
  deepSplit = 4, pamStage = FALSE, detectCutHeight = 0.995, minModuleSize = 50,
  mergeCutHeight = 0.2, saveConsensusTOMs = TRUE, ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  if (store_tom_in_seurat && !saveConsensusTOMs) {
    stop('store_tom_in_seurat = TRUE requires saveConsensusTOMs = TRUE.')
  }

  # suffix for the tom
  if (is.null(tom_name)) {
    tom_name <- gsub(' ', '_', wgcna_name)
  }

  # tom file name
  renamed <- paste0(tom_outdir, '/', tom_name, '_TOM.rda')
  if (file.exists(renamed)) {
    if (overwrite_tom) {
      warning(paste0('Overwriting TOM ', renamed))
    } else {
      stop(paste0('TOM ', renamed, ' already exists. Set overwrite_tom = TRUE or change tom_name to proceed.'))
    }
  }

  # constructing network on multiple datasets (consensus WGCNA)
  if (consensus) {
    multiExpr <- GetMultiExpr(seurat_obj, wgcna_name)
    checkSets(multiExpr) # check data size

    # constructing network on a single dataset
  } else {
    # get datExpr from seurat object
    datExpr <- GetDatExpr(seurat_obj, wgcna_name)

    if (is.null(wgcna_name)) {
      wgcna_name <- 'all'
    }

    nSets <- 1
    setLabels <- gsub(' ', '_', wgcna_name)
    shortLabels <- setLabels
    multiExpr <- list()
    multiExpr[[wgcna_name]] <- list(data = datExpr)
    checkSets(multiExpr) # check data size
  }

  # make output dir for the TOM
  if (!dir.exists(tom_outdir)) {
    dir.create(tom_outdir)
  }

  if (is.null(soft_power) && !consensus) {
    soft_power <- GetPowerTable(seurat_obj, wgcna_name) %>%
      subset(SFT.R.sq >= 0.8 & Power > min_power) %>%
      .$Power %>%
      min()
    cat(paste0('Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = ', soft_power, '\n'))
  } else if (is.null(soft_power)) {
    power_tables <- GetPowerTable(seurat_obj, wgcna_name) %>% dplyr::group_split(group)
    soft_power <- sapply(power_tables, function(power_table) {
      power_table %>%
        subset(SFT.R.sq >= 0.8 & Power > min_power) %>%
        .$Power %>%
        min()
    })
    cat(paste0('Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = c(', paste0(soft_power, collapse = ','), ')\n'))
  }

  # WGCNA writes temporary files (individual TOMs, disk cache) to the working directory by
  # default; keep them in tom_outdir with the tom_name prefix so that parallel runs from the
  # same directory don't overwrite each other's files. Users can still override them via ...
  wgcna_args <- list(...)
  if (!('individualTOMFileNames' %in% names(wgcna_args))) {
    wgcna_args$individualTOMFileNames <- file.path(
      tom_outdir, paste0(tom_name, '_individualTOM-Set%s-Block%b.RData')
    )
  }
  if (!('cacheDir' %in% names(wgcna_args))) {
    wgcna_args$cacheDir <- tom_outdir
  }
  if (!('cacheBase' %in% names(wgcna_args))) {
    wgcna_args$cacheBase <- paste0('.', tom_name, '_blockConsModsCache')
  }

  # construct the network
  net <- do.call(WGCNA::blockwiseConsensusModules, c(list(
    multiExpr,
    power = soft_power,
    blocks = blocks,
    maxBlockSize = maxBlockSize, ## This should be set to a smaller size if the user has limited RAM
    randomSeed = randomSeed,
    corType = corType,
    consensusQuantile = consensusQuantile,
    networkType = networkType,
    TOMType = TOMType,
    TOMDenom = TOMDenom,
    scaleTOMs = scaleTOMs,
    calibrationQuantile = calibrationQuantile,
    sampleForCalibration = sampleForCalibration,
    sampleForCalibrationFactor = sampleForCalibrationFactor,
    useDiskCache = useDiskCache, chunkSize = chunkSize,
    deepSplit = deepSplit,
    pamStage = pamStage,
    detectCutHeight = detectCutHeight, minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    saveConsensusTOMs = saveConsensusTOMs,
    consensusTOMFilePattern = paste0(tom_outdir, '/', tom_name, '_block.%b.rda')
  ), wgcna_args))

  # rename consensusTOM file:
  file.rename(paste0(tom_outdir, '/', tom_name, '_block.1.rda'), renamed)

  # add network parameters to the Seurat object:
  params <- list(
    power = soft_power,
    blocks = blocks,
    maxBlockSize = maxBlockSize, ## This should be set to a smaller size if the user has limited RAM
    randomSeed = randomSeed,
    corType = corType,
    consensusQuantile = consensusQuantile,
    networkType = networkType,
    TOMType = TOMType,
    TOMDenom = TOMDenom,
    scaleTOMs = scaleTOMs, calibrationQuantile = calibrationQuantile,
    sampleForCalibration = sampleForCalibration, sampleForCalibrationFactor = sampleForCalibrationFactor,
    useDiskCache = useDiskCache, chunkSize = chunkSize,
    deepSplit = deepSplit,
    pamStage = pamStage,
    detectCutHeight = detectCutHeight, minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    saveConsensusTOMs = saveConsensusTOMs
  )

  # add parameters:
  seurat_obj <- SetWGCNAParams(seurat_obj, params, wgcna_name)

  # append working directory to the TOM file so it has the full path:
  # absolute path to the TOM file (works for both relative and absolute tom_outdir)
  net$TOMFiles <- normalizePath(renamed, mustWork = FALSE)

  # add network to seurat obj
  # hdWGCNA identifies modules by colour; numericLabels = TRUE (passed via ...) makes WGCNA
  # return numbers instead, which breaks plotting. Convert them (0 -> 'grey').
  if (is.numeric(net$colors)) {
    message('Converting numeric module labels to colours (numericLabels is not supported).')
    gene_names <- names(net$colors)
    net$colors <- WGCNA::labels2colors(net$colors)
    names(net$colors) <- gene_names
  }

  seurat_obj <- SetNetworkData(seurat_obj, net, wgcna_name)

  # set the modules df in the Seurat object
  mods <- GetNetworkData(seurat_obj, wgcna_name)$colors
  seurat_obj <- SetModules(
    seurat_obj,
    modules = data.frame(
      'gene_name' = names(mods),
      'module' = factor(mods, levels = unique(mods)),
      'color' = mods
    ), wgcna_name
  )

  if (store_tom_in_seurat) {
    load(renamed)
    TOM_mat <- as.matrix(consTomDS)
    gene_names <- names(mods)
    rownames(TOM_mat) <- gene_names
    colnames(TOM_mat) <- gene_names
    seurat_obj <- SetTOM(seurat_obj, TOM_mat, wgcna_name)
  }

  seurat_obj
}


#' PlotDendrogram
#'
#' Plot WGCNA dendrogram
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
PlotDendrogram <- function(
  seurat_obj, groupLabels = 'Module colors', wgcna_name = NULL,
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = '', ...
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # get WGCNA network and module data
  net <- GetNetworkData(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)

  # plot dendrogram
  WGCNA::plotDendroAndColors(
    net$dendrograms[[1]],
    as.character(modules$color),
    groupLabels = groupLabels,
    dendroLabels = dendroLabels,
    hang = hang,
    addGuide = addGuide,
    guideHang = guideHang,
    main = main,
    ...
  )
}
