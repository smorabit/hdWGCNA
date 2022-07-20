#' ConstructMetacells
#'
#' This function takes a Seurat object and constructs averaged 'metacells' based
#' on neighboring cells.
#' @param seurat_obj A Seurat object
#' @param k Number of nearest neighbors to aggregate. Default = 50
#' @param name A string appended to resulting metalcells. Default = 'agg'
#' @param reduction A dimensionality reduction stored in the Seurat object. Default = 'umap'
#' @param assay Assay to extract data for aggregation. Default = 'RNA'
#' @param slot Slot to extract data for aggregation. Default = 'data'
#' @param return_metacell Logical to determine if we return the metacell seurat object (TRUE), or add it to the misc in the original Seurat object (FALSE). Default to FALSE.
#' @param mode determines how to make gene expression profiles for metacells from their constituent single cells. Options are "average" or "sum".
#' @param max_shared the maximum number of cells to be shared across two metacells
#' @param target_metacells the maximum target number of metacells to construct
#' @param max_iter the maximum number of iterations in the metacells bootstrapping loop
#' @param max_shared the maximum number of cells to be shared across two metacells
#' @param verbose logical indicating whether to print additional information
#' @param wgcna_name name of the WGCNA experiment
#' @keywords scRNA-seq
#' @export
#' @examples
#' ConstructMetacells
ConstructMetacells <- function(
  seurat_obj, name='agg', ident.group='seurat_clusters', k=25,
  reduction='umap', assay='RNA',
  cells.use = NULL, # if we don't want to use all the cells to make metacells, good for train/test split
  slot='counts',  meta=NULL, return_metacell=FALSE,
  mode = 'average', max_shared=15,
  target_metacells=1000,
  max_iter=5000,
  verbose=FALSE,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # check reduction
  if(!(reduction %in% names(seurat_obj@reductions))){
    stop(paste0("Invalid reduction (", reduction, "). Reductions in Seurat object: ", paste(names(seurat_obj@reductions), collapse=', ')))
  }

  # check assay
  if(!(assay %in% names(seurat_obj@assays))){
    stop(paste0("Invalid assay (", assay, "). Assays in Seurat object: ", paste(names(seurat_obj@assays), collapse=', ')))
  }

  # check slot
  if(!(slot %in% c('counts', 'data', 'scale.data'))){
    stop(paste0("Invalid slot (", slot, "). Valid options for slot: counts, data, scale.data "))
  }

  # subset seurat object by selected cells:
  if(!is.null(cells.use)){
    seurat_full <- seurat_obj
    seurat_obj <- seurat_obj[,cells.use]
  }


  reduced_coordinates <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
  nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), size = 1,
      replace = FALSE)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  it <- 0
  k2 <- k * 2
  get_shared <- function(other, this_choice) {
      k2 - length(union(cell_sample[other, ], this_choice))
  }
  while (length(good_choices) > 0 & length(chosen) < target_metacells & it < max_iter) {
      it <- it + 1
      choice <- sample(seq_len(length(good_choices)), size = 1,
          replace = FALSE)
      new_chosen <- c(chosen, good_choices[choice])
      good_choices <- good_choices[good_choices != good_choices[choice]]
      cell_sample <- nn_map[new_chosen, ]
      others <- seq_len(nrow(cell_sample) - 1)
      this_choice <- cell_sample[nrow(cell_sample), ]
      shared <- sapply(others, get_shared, this_choice = this_choice)

      # if (max(shared) < 0.9 * k) { # old way of doing it
      if(max(shared) <= max_shared){
          chosen <- new_chosen
      }
  }

  shared_old <- shared
  cell_sample <- nn_map[chosen, ]

  #
  # combs <- combn(nrow(cell_sample), 2)

  combs <- tryCatch(
    {combn(nrow(cell_sample), 2)},
    error = function(cond){return(NA)}
  )
  if(any(is.na(combs))){
    warning('Metacell failed')
    return(NULL)
  }


  shared <- apply(combs, 2, function(x) {
      k2 - length(unique(as.vector(cell_sample[x, ])))
  })

  if(verbose){
    message(paste0("Overlap QC metrics:\nCells per bin: ",
        k, "\nMaximum shared cells bin-bin: ", max(shared),
        "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ",
        median(shared)))
    if (mean(shared)/k > 0.1)
        warning("On average, more than 10% of cells are shared between paired bins.")
  }

  # get original expression matrix
  exprs_old <- GetAssayData(seurat_obj, assay=assay, slot=slot)

  # groups of cells to combine
  mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in%
      cell_sample[x, , drop = FALSE])
  # mask <- mask[,which(shared_old <= max_shared)]
  # cell_sample <- cell_sample[which(shared_old <= max_shared),]
  mask <- Matrix::Matrix(mask)

  # average or sum expression?
  new_exprs <- (exprs_old %*% mask)
  if(mode == 'average'){
    new_exprs <- new_exprs / k
  }
  colnames(new_exprs) <- paste0(name, '_', 1:ncol(new_exprs))
  rownames(cell_sample) <- paste0(name, '_', 1:ncol(new_exprs))
  colnames(cell_sample) <- paste0('knn_', 1:ncol(cell_sample))

  # make seurat obj:
  metacell_obj <- CreateSeuratObject(
    counts = new_exprs,
    assay = assay
  )

  # calculate stats:
  # shared <- shared[shared <= max_shared]
  max_shared <- max(shared)
  median_shared <- median(shared)
  mean_shared <- mean(shared)

  # calculate matrix density:
  new_exprs[new_exprs > 0] <- 1
  density <- sum(Matrix::colSums(new_exprs) / (nrow(new_exprs)*ncol(new_exprs)))
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
  if(!is.null(meta)){
    meta_names <- names(meta)
    for(x in meta_names){
      metacell_obj@meta.data[[x]] <- meta[[x]]
    }
  } else(
    warning('meta not found')
  )

  # add seurat metacell object to the main seurat object:
  if(return_metacell){
    out <- metacell_obj
  } else{

    # revert to full seurat object if we subsetted earlier
    if(!is.null(cells.use)){
      seurat_obj <- seurat_full
    }

    # add seurat metacell object to the main seurat object:
    seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)

    # add other info
    seurat_obj <- SetWGCNAParams(
      seurat_obj, params = list(
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

#' MetacellsByGroups
#'
#' This function takes a Seurat object and constructs averaged 'metacells' based
#' on neighboring cells in provided groupings, such as cluster or cell type.
#' @param seurat_obj A Seurat object
#' @param group.by A character vector of Seurat metadata column names representing groups for which metacells will be computed.
#' @param k Number of nearest neighbors to aggregate. Default = 50
#' @param name A string appended to resulting metalcells. Default = 'agg'
#' @param reduction A dimensionality reduction stored in the Seurat object. Default = 'pca'
#' @param assay Assay to extract data for aggregation. Default = 'RNA'
#' @param slot Slot to extract data for aggregation. Default = 'data'
#' @param mode determines how to make gene expression profiles for metacells from their constituent single cells. Options are "average" or "sum".
#' @param min_cells the minimum number of cells in a particular grouping to construct metacells
#' @param max_shared the maximum number of cells to be shared across two metacells
#' @param target_metacells the maximum target number of metacells to construct
#' @param max_iter the maximum number of iterations in the metacells bootstrapping loop
#' @param verbose logical indicating whether to print additional information
#' @param wgcna_name name of the WGCNA experiment
#' @keywords scRNA-seq
#' @export
#' @examples
#' MetacellsByGroups
MetacellsByGroups <- function(
  seurat_obj, group.by=c('seurat_clusters'),
  ident.group='seurat_clusters',
  k=25, reduction='pca', assay=NULL,
  cells.use = NULL, # if we don't want to use all the cells to make metacells, good for train/test split
  slot='counts', mode = 'average', min_cells=100,
  max_shared=15,
  target_metacells=1000,
  max_iter=5000, verbose=FALSE, wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # check group.by for invalid characters:
  if(any(grepl('#', group.by))){
    stop('Invalid character # found in group.by, please re-name the group.')
  }

  if(!(ident.group %in% group.by)){
    stop('ident.group must be in group.by')
  }

  if(!(mode %in% c('sum', 'average'))){
    stop('Invalid choice for mode. Mode can be either sum or average.')
  }

  # check reduction
  if(!(reduction %in% names(seurat_obj@reductions))){
    stop(paste0("Invalid reduction (", reduction, "). Reductions in Seurat object: ", paste(names(seurat_obj@reductions), collapse=', ')))
  }

  # check assay:
  if(is.null(assay)){
    assay <- DefaultAssay(seurat_obj)
  } else if(!(assay %in% names(seurat_obj@assays))){
    stop(paste0('Assay ', assay, ' not found in seurat_obj. Select a valid assay: ', paste0(names(seurat_obj@assays), collapse = ', ')))
  }

  # subset seurat object by seleted cells:
  if(!is.null(cells.use)){
    seurat_full <- seurat_obj
    seurat_obj <- seurat_obj[,cells.use]
  }

  # setup grouping variables
  if(length(group.by) > 1){
    seurat_meta <- seurat_obj@meta.data[,group.by]
    for(col in colnames(seurat_meta)){
      seurat_meta[[col]] <- as.character(seurat_meta[[col]])
    }
    seurat_obj$metacell_grouping <- apply(seurat_meta, 1, paste, collapse='#')
  } else {
    seurat_obj$metacell_grouping <- as.character(seurat_obj@meta.data[[group.by]])
  }
  groupings <- unique(seurat_obj$metacell_grouping)
  groupings <- groupings[order(groupings)]

  # remove groups that are too small:
  group_counts <- table(seurat_obj$metacell_grouping) >= min_cells
  warning(paste0("Removing the following groups that did not meet min_cells: ", paste(names(group_counts)[group_counts], collapse=', ')))
  groupings <- groupings[table(seurat_obj$metacell_grouping) >= min_cells]

  if(length(groupings) == 0 ){
    stop("No groups met the min_cells requirement.")
  }

  # unique meta-data for each group
  meta_df <- as.data.frame(do.call(rbind, strsplit(groupings, '#')))
  colnames(meta_df) <- group.by

  # list of meta-data to pass to each metacell seurat object
  meta_list <- lapply(1:nrow(meta_df), function(i){
    x <- list(as.character(meta_df[i,]))[[1]]
    names(x) <- colnames(meta_df)
    x
  })

  # split seurat obj by groupings
  seurat_list <- lapply(groupings, function(x){seurat_obj[,seurat_obj$metacell_grouping == x]})
  names(seurat_list) <- groupings

  # construct metacells
  metacell_list <- mapply(
    ConstructMetacells,
    seurat_obj = seurat_list,
    name = groupings,
    meta = meta_list,
    MoreArgs = list(k=k, reduction=reduction, assay=assay, slot=slot, return_metacell=TRUE, mode=mode, max_shared=max_shared, max_iter=max_iter, target_metacells=target_metacells, verbose=verbose, wgcna_name=wgcna_name)
  )
  names(metacell_list) <- groupings


  # remove NULL
  remove <- which(sapply(metacell_list, is.null))
  if(length(remove) > 1){
    metacell_list <- metacell_list[-remove]
  }
  # get the run stats:
  run_stats <- as.data.frame(do.call(rbind, lapply(metacell_list, function(x){x@misc$run_stats})))
  rownames(run_stats) <- 1:nrow(run_stats)
  for(i in 1:length(group.by)){
    run_stats[[group.by[i]]] <- do.call(rbind, strsplit(run_stats$name, '#'))[,i]
  }

  # combine metacell objects
  metacell_obj <- merge(metacell_list[[1]], metacell_list[2:length(metacell_list)])

  # set idents for metacell object:
  Idents(metacell_obj) <- metacell_obj@meta.data[[ident.group]]

  # revert to full seurat object if we subsetted earlier
  if(!is.null(cells.use)){
    seurat_obj <- seurat_full
  }

  # add seurat metacell object to the main seurat object:
  seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)

  # add other info
  seurat_obj <- SetWGCNAParams(
    seurat_obj, params = list(
      'metacell_k' = k,
      'metacell_reduction' = reduction,
      'metacell_slot' = slot,
      'metacell_assay' = assay,
      'metacell_stats' = run_stats
    ),
    wgcna_name
  )

  seurat_obj
}
