#' MetacellsByGroups
#'
#' This function takes a Seurat object and constructs 'metacells' by aggregating
#' neighboring cells within the provided groupings, such as cell type and sample.
#'
#' @return The input Seurat object with the metacell Seurat object and metacell run statistics stored in \code{seurat_obj@misc[[wgcna_name]]}; retrieve the metacells with GetMetacellObject().
#'
#' @param seurat_obj A Seurat object.
#' @param group.by Required. A character vector of Seurat metadata column names; metacells are built separately within each combination of these columns. Typically the cell population and the biological replicate, e.g. \code{c('cell_type', 'Sample')}, so that no metacell mixes cells from different samples.
#' @param ident.group A column name in group.by used to set the Idents of the metacell object. Default = NULL to use the first column of group.by.
#' @param k Number of nearest neighbors aggregated into each metacell. Default = 50. Use a smaller k (e.g. 25) for cell populations with few cells per group.
#' @param reduction Required. The name of the dimensionality reduction used to find nearest neighbors, typically the integrated embedding (e.g. 'harmony').
#' @param dims A vector represnting the dimensions of the reduction to use. Either specify the names of the dimensions or the indices. Default = NULL to include all dims.
#' @param assay Assay to extract data for aggregation. Default = NULL to use the default assay.
#' @param slot Slot to extract data for aggregation. Default = 'counts'. Slot is used with Seurat v4 instead of layer.
#' @param layer Layer to extract data for aggregation. Default = 'counts'. Layer is used with Seurat v5 instead of slot.
#' @param mode How to combine the expression of the cells in a metacell: 'sum' (default) or 'average'.
#' @param cells.use A character vector of cell barcodes to use for metacell construction. Default = NULL to use all cells.
#' @param min_cells The minimum number of cells in a group; smaller groups are skipped. Default = 100
#' @param max_shared The maximum number of cells shared between any two metacells. Default = NULL to use k / 5 (10 for k = 50, 5 for k = 25).
#' @param target_metacells The maximum number of metacells to construct per group. Default = 100
#' @param max_iter the maximum number of iterations in the metacells bootstrapping loop
#' @param verbose logical indicating whether to print additional information
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#'
#' @details
#' MetacellsByGroups merges transcriptomically similar cells into "metacells".
#' Given a dimensionally-reduced representation of the input dataset, this algorithm
#' first uses KNN to identify similar cells. A bootstrapped sampling procedure is then
#' used to group together similar cells until convergence is reached. Importantly,
#' this procedure is done in a context-specific manner based on the provided group.by parameters.
#' Typically this means that metacells will be constructed separately for each biological
#' replicate, cell type or cell state, disease condition, etc. The metacell representation is
#' considerably less sparse than the original single-cell dataset, which is preferable for
#' co-expression network analysis or othter analyses that rely on correlations.
#'
#' @export
MetacellsByGroups <- function(
  seurat_obj,
  group.by,
  ident.group = NULL,
  k = 50,
  reduction,
  dims = NULL,
  assay = NULL,
  slot = 'counts',
  layer = 'counts',
  mode = 'sum',
  cells.use = NULL,
  min_cells = 100,
  max_shared = NULL,
  target_metacells = 100,
  max_iter = 5000,
  verbose = FALSE,
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # group.by and reduction have no defaults: the right choice depends on the dataset
  if (missing(group.by)) {
    stop(paste0(
      'group.by is required: the metadata columns to build metacells within, typically ',
      "the cell population and the sample, e.g. group.by = c('cell_type', 'Sample')."
    ))
  }
  if (missing(reduction)) {
    stop(paste0(
      'reduction is required: the embedding used to find nearest neighbors, typically ',
      "the integrated one, e.g. reduction = 'harmony'. Reductions in this object: ",
      paste(names(seurat_obj@reductions), collapse = ', ')
    ))
  }
  missing_cols <- setdiff(group.by, colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(paste0(
      'group.by column(s) not found in the metadata: ',
      paste(missing_cols, collapse = ', ')
    ))
  }

  # check group.by for invalid characters:
  if (any(grepl('#', group.by))) {
    stop('Invalid character # found in group.by, please re-name the group.')
  }

  # check ident.group
  if (is.null(ident.group)) {
    ident.group <- group.by[1]
  }
  if (!(ident.group %in% group.by)) {
    stop('ident.group must be in group.by')
  }

  # check mode
  if (!(mode %in% c('sum', 'average'))) {
    stop('Invalid choice for mode. Mode can be either sum or average.')
  }

  # check reduction
  if (!(reduction %in% names(seurat_obj@reductions))) {
    stop(paste0('Invalid reduction (', reduction, '). Reductions in Seurat object: ', paste(names(seurat_obj@reductions), collapse = ', ')))
  }

  # check assay:
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  } else if (!(assay %in% names(seurat_obj@assays))) {
    stop(paste0('Assay ', assay, ' not found in seurat_obj. Select a valid assay: ', paste0(names(seurat_obj@assays), collapse = ', ')))
  }

  # check slot/layer:
  if (CheckSeurat5()) {
    if (!(layer %in% c('counts', 'data', 'scale.data'))) {
      stop('Invalid input for layer. Valid choices are counts, data, scale.data.')
    } else {
      # check the shape of the matrix
      layer_dim <- dim(SeuratObject::LayerData(seurat_obj, assay = assay, layer = layer))
      if (any(layer_dim == 0)) {
        stop(paste(c('Selected layer ', layer, ' not found in this assay.')))
      }
    }
  } else {
    if (!(slot %in% c('counts', 'data', 'scale.data'))) {
      stop('Invalid input for slot. Valid choices are counts, data, scale.data.')
    } else {
      # check the shape of the matrix
      slot_dim <- dim(Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot))
      if (any(slot_dim == 0)) {
        stop(paste(c('Selected slot ', slot, ' not found in this assay.')))
      }
    }
  }

  # check that k > min_cells
  if (min_cells < k) {
    warning('min_cells is smaller than k, this may result in downstream errors if very small groups are allowed.')
  }

  # check that max_shared is valid:
  if (is.null(max_shared)) {
    max_shared <- default_max_shared(k)
  }
  if (max_shared < 0) {
    max_shared <- 0
    warning(paste0('max_shared specified (', max_shared, ') is too low, setting max_shared <- 0'))
  }

  # subset seurat object by seleted cells:
  if (!is.null(cells.use)) {
    seurat_full <- seurat_obj
    seurat_obj <- seurat_obj[, cells.use]
  }

  # setup grouping variables
  if (length(group.by) > 1) {
    seurat_meta <- seurat_obj@meta.data[, group.by]
    for (col in colnames(seurat_meta)) {
      seurat_meta[[col]] <- as.character(seurat_meta[[col]])
    }
    seurat_obj$metacell_grouping <- apply(seurat_meta, 1, paste, collapse = '#')
  } else {
    seurat_obj$metacell_grouping <- as.character(seurat_obj@meta.data[[group.by]])
  }
  groupings <- unique(seurat_obj$metacell_grouping)
  groupings <- groupings[order(groupings)]

  # remove groups that are too small:
  group_counts <- table(seurat_obj$metacell_grouping) < min_cells
  if (any(group_counts)) {
    warning(paste0('Removing the following groups that did not meet min_cells: ', paste(names(group_counts)[group_counts], collapse = ', ')))
  }
  group_sizes <- table(seurat_obj$metacell_grouping)
  groupings <- groupings[groupings %in% names(group_sizes)[group_sizes >= min_cells]]

  if (length(groupings) == 0) {
    stop('No groups met the min_cells requirement.')
  }

  # unique meta-data for each group
  meta_df <- as.data.frame(do.call(rbind, strsplit(groupings, '#')))
  colnames(meta_df) <- group.by

  # list of meta-data to pass to each metacell seurat object
  meta_list <- lapply(seq_len(nrow(meta_df)), function(i) {
    x <- list(as.character(meta_df[i, ]))[[1]]
    names(x) <- colnames(meta_df)
    x
  })

  # split seurat obj by groupings
  # split a lightweight copy: graphs and neighbors are not needed here, and subsetting them
  # fails for graphs without dimnames (e.g. objects from sketch / ProjectData workflows)
  split_obj <- seurat_obj
  split_obj@graphs <- list()
  split_obj@neighbors <- list()
  seurat_list <- lapply(groupings, function(x) {
    split_obj[, split_obj$metacell_grouping == x]
  })
  names(seurat_list) <- groupings

  # construct metacells
  metacell_list <- mapply(
    ConstructMetacells,
    seurat_obj = seurat_list,
    name = groupings,
    meta = meta_list,
    MoreArgs = list(
      k = k,
      reduction = reduction,
      dims = dims,
      assay = assay,
      slot = slot,
      layer = layer,
      return_metacell = TRUE,
      mode = mode,
      max_shared = max_shared,
      max_iter = max_iter,
      target_metacells = target_metacells,
      verbose = verbose,
      wgcna_name = wgcna_name
    )
  )
  names(metacell_list) <- groupings

  # remove groups where metacell construction failed
  remove <- which(vapply(metacell_list, is.null, logical(1)))
  if (length(remove) == length(metacell_list)) {
    stop(paste0(
      'Metacell construction failed for all groups. Each group needs enough cells to form ',
      'metacells: try lowering k, raising max_shared, or making sure min_cells is well above k.'
    ))
  }
  if (length(remove) >= 1) {
    warning(paste0(
      'Metacell construction failed for ', length(remove), ' group(s), which were removed: ',
      paste0(names(metacell_list)[remove], collapse = ', '),
      '. Try lowering k, raising max_shared, or raising min_cells to exclude small groups.'
    ))
    metacell_list <- metacell_list[-remove]
  }

  # get the run stats:
  run_stats <- as.data.frame(do.call(rbind, lapply(metacell_list, function(x) {
    x@misc$run_stats
  })))
  rownames(run_stats) <- seq_len(nrow(run_stats))
  for (i in seq_along(group.by)) {
    run_stats[[group.by[i]]] <- do.call(rbind, strsplit(as.character(run_stats$name), '#'))[, i]
  }
  message(metacell_summary(run_stats, k = k))

  # combine metacell objects
  if (length(metacell_list) > 1) {
    metacell_obj <- merge(metacell_list[[1]], metacell_list[2:length(metacell_list)])

    # join the layers of Seurat v5 assays (v3 assays have no layers to join)
    if (inherits(metacell_obj[[assay]], 'Assay5')) {
      metacell_obj <- SeuratObject::JoinLayers(metacell_obj)
    }
  } else {
    metacell_obj <- metacell_list[[1]]
  }

  # set idents for metacell object:
  Idents(metacell_obj) <- metacell_obj@meta.data[[ident.group]]

  # revert to full seurat object if we subsetted earlier
  if (!is.null(cells.use)) {
    seurat_obj <- seurat_full
  }

  # add seurat metacell object to the main seurat object:
  seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)

  # add other info
  param_list <- list(
    'metacell_k' = k,
    'metacell_reduction' = reduction,
    'metacell_assay' = assay,
    'metacell_stats' = run_stats
  )

  if (CheckSeurat5()) {
    param_list['metacell_layer'] <- layer
  } else {
    param_list['metacell_slot'] <- slot
  }

  seurat_obj <- SetWGCNAParams(
    seurat_obj,
    params = param_list, wgcna_name
  )

  # return the updated seurat object
  seurat_obj
}


#' ConstructMetacells
#'
#' This function takes a Seurat object and constructs averaged 'metacells' based
#' on neighboring cells.
#' @param seurat_obj A Seurat object.
#' @param name A string prefix for the resulting metacell names. Default = 'agg'
#' @param ident.group Not currently used by this function; kept for compatibility with MetacellsByGroups.
#' @param k Number of nearest neighbors aggregated into each metacell. Default = 50
#' @param reduction Required. The name of the dimensionality reduction used to find nearest neighbors.
#' @param dims A vector represnting the dimensions of the reduction to use. Either specify the names of the dimensions or the indices. Default = NULL to include all dims.
#' @param assay Assay to extract data for aggregation. Default = 'RNA'
#' @param cells.use A character vector of cell barcodes to use for metacell construction. Default = NULL to use all cells.
#' @param slot Slot to extract data for aggregation. Default = 'counts'. Slot is used with Seurat v4 instead of layer.
#' @param layer Layer to extract data for aggregation. Default = 'counts'. Layer is used with Seurat v5 instead of slot.
#' @param meta A named list of metadata values to add to every metacell (e.g. the grouping values). Default = NULL
#' @param return_metacell Logical to determine if we return the metacell seurat object (TRUE), or add it to the misc in the original Seurat object (FALSE). Default to FALSE.
#' @param mode How to combine the expression of the cells in a metacell: 'sum' (default) or 'average'.
#' @param max_shared The maximum number of cells shared between any two metacells. Default = NULL to use k / 5.
#' @param target_metacells The maximum number of metacells to construct. Default = 100
#' @param max_iter the maximum number of iterations in the metacells bootstrapping loop
#' @param verbose logical indicating whether to print additional information
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return If \code{return_metacell = TRUE}, the metacell Seurat object; otherwise the input Seurat object with the metacell object stored in \code{seurat_obj@misc[[wgcna_name]]} (retrieve with GetMetacellObject()). Returns NULL if metacell construction fails.
#' @keywords internal
#' @noRd
ConstructMetacells <- function(
  seurat_obj, name = 'agg', ident.group = 'seurat_clusters', k = 50,
  reduction,
  dims = NULL,
  assay = 'RNA',
  cells.use = NULL, # if we don't want to use all the cells to make metacells, good for train/test split
  slot = 'counts',
  layer = 'counts',
  meta = NULL, return_metacell = FALSE,
  mode = 'sum', max_shared = NULL,
  target_metacells = 100,
  max_iter = 5000,
  verbose = FALSE,
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  if (missing(reduction)) {
    stop(paste0(
      "reduction is required, e.g. reduction = 'harmony'. Reductions in this object: ",
      paste(names(seurat_obj@reductions), collapse = ', ')
    ))
  }
  if (is.null(max_shared)) {
    max_shared <- default_max_shared(k)
  }

  # check reduction
  if (!(reduction %in% names(seurat_obj@reductions))) {
    stop(paste0('Invalid reduction (', reduction, '). Reductions in Seurat object: ', paste(names(seurat_obj@reductions), collapse = ', ')))
  }

  # check assay
  if (!(assay %in% names(seurat_obj@assays))) {
    stop(paste0('Invalid assay (', assay, '). Assays in Seurat object: ', paste(names(seurat_obj@assays), collapse = ', ')))
  }

  # check slot
  if (!(slot %in% c('counts', 'data', 'scale.data'))) {
    stop(paste0('Invalid slot (', slot, '). Valid options for slot: counts, data, scale.data '))
  }

  # subset seurat object by selected cells:
  if (!is.null(cells.use)) {
    seurat_full <- seurat_obj
    seurat_obj <- seurat_obj[, cells.use]
  }

  # get the dim reduction
  check_reduction(seurat_obj, reduction)
  reduced_coordinates <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
  reduc_orig <- reduced_coordinates

  # subset the dims?
  if (!is.null(dims)) {
    # are these indices?
    if (is.numeric(dims)) {
      if (!all(dims %in% seq_len(ncol(reduced_coordinates)))) {
        stop('Invalid selection for dims. Check that the selected dims match the columns of the selected reduction.')
      }
    } else {
      if (!all(dims %in% colnames(reduced_coordinates))) {
        stop('Invalid selection for dims. Check that the selected dims match the column names of the selected reduction.')
      }
    }
    reduced_coordinates <- reduced_coordinates[, dims]
  }

  # run KNN on the chosen dim reduction
  nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))

  # set up variables for loop
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_along(good_choices),
    size = 1,
    replace = FALSE
  )
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  it <- 0
  k2 <- k * 2
  get_shared <- function(other, this_choice) {
    k2 - length(union(cell_sample[other, ], this_choice))
  }

  # loop to create metacells until convergencce
  while (length(good_choices) > 0 && length(chosen) < target_metacells && it < max_iter) {
    it <- it + 1
    choice <- sample(seq_along(good_choices),
      size = 1,
      replace = FALSE
    )
    new_chosen <- c(chosen, good_choices[choice])
    good_choices <- good_choices[good_choices != good_choices[choice]]
    cell_sample <- nn_map[new_chosen, ]
    others <- seq_len(nrow(cell_sample) - 1)
    this_choice <- cell_sample[nrow(cell_sample), ]
    shared <- sapply(others, get_shared, this_choice = this_choice)

    if (max(shared) <= max_shared) {
      chosen <- new_chosen
    }
  }

  shared_old <- shared
  cell_sample <- nn_map[chosen, ]

  if (length(chosen) <= 1) {
    warning('Metacell failed')
    return(NULL)
  }

  # get a list of the cell barcodes that have been merged:
  cells_merged <- apply(cell_sample, 1, function(x) {
    paste0(colnames(seurat_obj)[x], collapse = ',')
  })

  combs <- tryCatch(
    {
      combn(nrow(cell_sample), 2)
    },
    error = function(cond) {
      return(NA)
    }
  )
  if (any(is.na(combs))) {
    warning('Metacell failed')
    return(NULL)
  }

  shared <- apply(combs, 2, function(x) {
    k2 - length(unique(as.vector(cell_sample[x, ])))
  })

  if (verbose) {
    message(paste0(
      'Overlap QC metrics:\nCells per bin: ',
      k, '\nMaximum shared cells bin-bin: ', max(shared),
      '\nMean shared cells bin-bin: ', mean(shared), '\nMedian shared cells bin-bin: ',
      median(shared)
    ))
    if (mean(shared) / k > 0.1) {
      warning('On average, more than 10% of cells are shared between paired bins.')
    }
  }

  # get original expression matrix
  if (CheckSeurat5()) {
    exprs_old <- SeuratObject::LayerData(seurat_obj, assay = assay, layer = layer)
  } else {
    exprs_old <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  }

  # groups of cells to combine
  mask <- sapply(seq_len(nrow(cell_sample)), function(x) {
    seq_len(ncol(exprs_old)) %in%
      cell_sample[x, , drop = FALSE]
  })
  mask <- Matrix::Matrix(mask)

  # average or sum expression?
  new_exprs <- (exprs_old %*% mask)
  if (mode == 'average') {
    new_exprs <- new_exprs / k
  }
  colnames(new_exprs) <- paste0(name, '_', seq_len(ncol(new_exprs)))
  rownames(cell_sample) <- paste0(name, '_', seq_len(ncol(new_exprs)))
  colnames(cell_sample) <- paste0('knn_', seq_len(ncol(cell_sample)))

  # make seurat obj:
  metacell_obj <- CreateSeuratObject(
    counts = new_exprs,
    assay = assay
  )
  if (slot == 'scale.data') {
    metacell_obj <- SeuratObject::SetAssayData(
      metacell_obj,
      slot = slot,
      assay = assay,
      new.data = as.matrix(new_exprs)
    )
  }

  # add the cells merged info to the metacell obj
  metacell_obj$cells_merged <- as.character(cells_merged)

  # calculate stats:
  max_shared <- max(shared)
  median_shared <- median(shared)
  mean_shared <- mean(shared)

  # calculate matrix density:
  new_exprs[new_exprs > 0] <- 1
  density <- sum(Matrix::colSums(new_exprs) / (nrow(new_exprs) * ncol(new_exprs)))
  run_stats <- data.frame(
    name = name,
    max_shared = max_shared,
    mean_shared = mean_shared,
    median_shared = median_shared,
    density = density,
    n = ncol(new_exprs)
  )

  # add to metacell seurat obj
  metacell_obj@misc$run_stats <- run_stats

  # add meta-data:
  if (!is.null(meta)) {
    meta_names <- names(meta)
    for (x in meta_names) {
      metacell_obj@meta.data[[x]] <- meta[[x]]
    }
  } else {
    (
      warning('meta not found')
    )
  }

  # add seurat metacell object to the main seurat object:
  if (return_metacell) {
    out <- metacell_obj
  } else {
    # revert to full seurat object if we subsetted earlier
    if (!is.null(cells.use)) {
      seurat_obj <- seurat_full
    }

    # add seurat metacell object to the main seurat object:
    seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)

    # add other info
    seurat_obj <- SetWGCNAParams(
      seurat_obj,
      params = list(
        'metacell_k' = k,
        'metacell_reduction' = reduction,
        'metacell_slot' = slot,
        'metacell_assay' = assay
      ),
      wgcna_name
    )

    out <- seurat_obj
  }
  out
}


#' NormalizeMetacells
#'
#' Wrapper function to run Seurat's NormalizeData function on the metacell object.
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param ... Additional arguments passed to Seurat::NormalizeData
#' @return The input Seurat object with the normalized metacell object stored in \code{seurat_obj@misc[[wgcna_name]]}; retrieve with GetMetacellObject().
#' @export
NormalizeMetacells <- function(seurat_obj, wgcna_name = NULL, ...) {
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::NormalizeData(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}


#' ScaleMetacells
#'
#' Wrapper function to run Seurat's ScaleData function on the metacell object.
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param ... Additional arguments passed to Seurat::ScaleData
#' @return The input Seurat object with the scaled metacell object stored in \code{seurat_obj@misc[[wgcna_name]]}; retrieve with GetMetacellObject().
#' @export
ScaleMetacells <- function(seurat_obj, wgcna_name = NULL, ...) {
  if (!exists('features')) {
    features <- VariableFeatures(seurat_obj)
  }
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::ScaleData(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}


#' RunPCAMetacells
#'
#' Wrapper function to run Seurat's RunPCA function on the metacell object.
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @param ... Additional arguments passed to Seurat::RunPCA
#' @return The input Seurat object with the PCA-reduced metacell object stored in \code{seurat_obj@misc[[wgcna_name]]}; retrieve with GetMetacellObject().
#' @export
RunPCAMetacells <- function(seurat_obj, wgcna_name = NULL, ...) {
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::RunPCA(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}


#' RunHarmonyMetacells
#'
#' Wrapper function to run harmony's RunHarmony function on the metacell object.
#'
#' @param seurat_obj A Seurat object.
#' @param ... Additional arguments passed to harmony::RunHarmony
#' @return The input Seurat object with the Harmony-corrected metacell object stored in the active hdWGCNA experiment; retrieve with GetMetacellObject().
#' @export
RunHarmonyMetacells <- function(seurat_obj, ...) {
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- harmony::RunHarmony(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
  seurat_obj
}


#' RunUMAPMetacells
#'
#' Wrapper function to run Seurat's RunUMAP function on the metacell object.
#'
#' @param seurat_obj A Seurat object.
#' @param ... Additional arguments passed to Seurat::RunUMAP
#' @return The input Seurat object with the UMAP-reduced metacell object stored in the active hdWGCNA experiment; retrieve with GetMetacellObject().
#' @keywords scRNA-seq
#' @export
RunUMAPMetacells <- function(seurat_obj, ...) {
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::RunUMAP(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
  seurat_obj
}


#' DimPlotMetacells
#'
#' Wrapper function to run Seurat's DimPlot function on the metacell object.
#'
#' @param seurat_obj A Seurat object.
#' @param ... Additional arguments passed to Seurat::DimPlot
#' @return A ggplot object from Seurat::DimPlot of the metacell object in the active hdWGCNA experiment.
#' @keywords scRNA-seq
#' @export
DimPlotMetacells <- function(seurat_obj, ...) {
  Seurat::DimPlot(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
}


# Default overlap allowed between two metacells: k / 5 (10 for k = 50, 5 for k = 25).
default_max_shared <- function(k) {
  floor(k / 5)
}

# One-line report of how many metacells each group produced, naming groups that
# produced few, so it is clear when a population needs a smaller k.
metacell_summary <- function(run_stats, k, few = 10) {
  n <- run_stats$n
  out <- paste0(
    'Constructed ', sum(n), ' metacells (k = ', k, ') in ', length(n), ' group(s); ',
    'per group: median ', stats::median(n), ', range ', min(n), '-', max(n), '.'
  )
  small <- run_stats$name[n < few]
  if (length(small) > 0) {
    shown <- utils::head(small, 10)
    more <- if (length(small) > 10) paste0(' and ', length(small) - 10, ' more') else ''
    out <- paste0(
      out, '\n', length(small), ' group(s) produced fewer than ', few,
      ' metacells; consider a smaller k for these populations: ',
      paste(gsub('#', ' / ', shown), collapse = ', '), more, '.'
    )
  }
  out
}
