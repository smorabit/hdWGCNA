#' ConstructMetacells
#'
#' This function takes a Seurat object and constructs averaged 'metacells' based
#' on neighboring cells.
#' @param seurat_obj A Seurat object
#' @param k Number of nearest neighbors to aggregate. Default = 50
#' @param name A string appended to resulting metalcells. Default = 'agg'
#' @param reduction A dimensionality reduction stored in the Seurat object. Default = 'umap'
#' @param dims A vector represnting the dimensions of the reduction to use. Either specify the names of the dimensions or the indices. Default = NULL to include all dims.
#' @param assay Assay to extract data for aggregation. Default = 'RNA'
#' @param slot Slot to extract data for aggregation. Default = 'counts'. Slot is used with Seurat v4 instead of layer.
#' @param layer Layer to extract data for aggregation. Default = 'counts'. Layer is used with Seurat v5 instead of slot.
#' @param return_metacell Logical to determine if we return the metacell seurat object (TRUE), or add it to the misc in the original Seurat object (FALSE). Default to FALSE.
#' @param mode determines how to make gene expression profiles for metacells from their constituent single cells. Options are "average" or "sum".
#' @param max_shared the maximum number of cells to be shared across two metacells
#' @param target_metacells the maximum target number of metacells to construct
#' @param max_iter the maximum number of iterations in the metacells bootstrapping loop
#' @param max_shared the maximum number of cells to be shared across two metacells
#' @param verbose logical indicating whether to print additional information
#' @param wgcna_name name of the WGCNA experiment
#' @export
#' ConstructMetacells
ConstructMetacells <- function(
  seurat_obj, name='agg', ident.group='seurat_clusters', k=25,
  reduction='pca', 
  dims = NULL,
  assay='RNA',
  cells.use = NULL, # if we don't want to use all the cells to make metacells, good for train/test split
  slot='counts',  
  layer='counts',
  meta=NULL, return_metacell=FALSE,
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

  # get the dim reduction
  reduced_coordinates <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
  reduc_orig <- reduced_coordinates

  # subset the dims?
  if(!is.null(dims)){
    # are these indices?
    if(is.numeric(dims)){
      if(!all(dims %in% 1:ncol(reduced_coordinates))){
        stop("Invalid selection for dims. Check that the selected dims match the columns of the selected reduction.")
      }
    } else{
      if(!all(dims %in% colnames(reduced_coordinates))){
        stop("Invalid selection for dims. Check that the selected dims match the column names of the selected reduction.")
      }
    }
    reduced_coordinates <- reduced_coordinates[,dims]
  }

  # run KNN on the chosen dim reduction
  nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  
  # set up variables for loop
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

  # loop to create metacells until convergencce
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

      if(max(shared) <= max_shared){
          chosen <- new_chosen
      }
  }

  shared_old <- shared
  cell_sample <- nn_map[chosen, ]

   if(length(chosen) <= 1){
      warning('Metacell failed')
      return(NULL)
    }

  # get a list of the cell barcodes that have been merged:
  cells_merged <- apply(cell_sample, 1, function(x){
    paste0(colnames(seurat_obj)[x], collapse=',')
  })

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
  if(CheckSeurat5()){
    exprs_old <- SeuratObject::LayerData(seurat_obj, assay=assay, layer=layer)
  } else{
    exprs_old <- Seurat::GetAssayData(seurat_obj, assay=assay, slot=slot)
  }

  # groups of cells to combine
  mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in%
      cell_sample[x, , drop = FALSE])
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
  if(slot == 'scale.data'){
    metacell_obj <- SeuratObject::SetAssayData(
      metacell_obj,
      slot=slot,
      assay=assay,
      new.data=as.matrix(new_exprs)
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
#' 
#' @return seurat_obj with a metacell seurat_obj stored in the specified WGCNA experiment
#' 
#' @param seurat_obj A Seurat object
#' @param group.by A character vector of Seurat metadata column names representing groups for which metacells will be computed.
#' @param k Number of nearest neighbors to aggregate. Default = 50
#' @param name A string appended to resulting metalcells. Default = 'agg'
#' @param reduction A dimensionality reduction stored in the Seurat object. Default = 'pca'
#' @param dims A vector represnting the dimensions of the reduction to use. Either specify the names of the dimensions or the indices. Default = NULL to include all dims.
#' @param assay Assay to extract data for aggregation. Default = 'RNA'
#' @param slot Slot to extract data for aggregation. Default = 'counts'. Slot is used with Seurat v4 instead of layer.
#' @param layer Layer to extract data for aggregation. Default = 'counts'. Layer is used with Seurat v5 instead of slot.
#' @param mode determines how to make gene expression profiles for metacells from their constituent single cells. Options are "average" or "sum".
#' @param min_cells the minimum number of cells in a particular grouping to construct metacells
#' @param max_shared the maximum number of cells to be shared across two metacells
#' @param target_metacells the maximum target number of metacells to construct
#' @param max_iter the maximum number of iterations in the metacells bootstrapping loop
#' @param verbose logical indicating whether to print additional information
#' @param wgcna_name name of the WGCNA experiment
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
#' @import Seurat
#' @export
MetacellsByGroups <- function(
  seurat_obj, group.by=c('seurat_clusters'),
  ident.group='seurat_clusters',
  k=25, 
  reduction='pca', 
  dims=NULL,
  assay=NULL,
  slot='counts', 
  layer='counts',
  mode = 'average', 
  cells.use = NULL, 
  min_cells=100,
  max_shared=15,
  target_metacells=1000,
  max_iter=5000, 
  verbose=FALSE, 
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # check group.by for invalid characters:
  if(any(grepl('#', group.by))){
    stop('Invalid character # found in group.by, please re-name the group.')
  }

  # check ident.group
  if(!(ident.group %in% group.by)){
    stop('ident.group must be in group.by')
  }

  # check mode
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

  # check slot/layer:
  if(!(slot %in% c('counts', 'data', 'scale.data'))){
    stop('Invalid input for slot. Valid choices are counts, data, scale.data.')
  } else{

    # check the shape of the slot
    slot_dim <- dim(Seurat::GetAssayData(seurat_obj, assay=assay, slot=slot))
    if(any(slot_dim) == 0){
      stop(paste(c("Selected slot ", slot, " not found in this assay.")))
    }
  }

  # check that k > min_cells 
  if(min_cells < k ){
    warning("min_cells is smaller than k, this may result in downstream errors if very small groups are allowed.")
  }

  # check that max_shared is valid:
  if(max_shared < 0){
    max_shared <- 0
    warning(paste0("max_shared specified (", max_shared, ') is too low, setting max_shared <- 0'))
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
  group_counts <- table(seurat_obj$metacell_grouping) < min_cells
  if(any(group_counts)){
    warning(paste0("Removing the following groups that did not meet min_cells: ", paste(names(group_counts)[group_counts], collapse=', ')))
  }
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
    MoreArgs = list(
      k=k, 
      reduction=reduction, 
      dims=dims,
      assay=assay, 
      slot=slot, 
      layer=layer,
      return_metacell=TRUE, 
      mode=mode,
      max_shared=max_shared, 
      max_iter=max_iter, 
      target_metacells=target_metacells, 
      verbose=verbose, 
      wgcna_name=wgcna_name
    )
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
    run_stats[[group.by[i]]] <- do.call(rbind, strsplit(as.character(run_stats$name), '#'))[,i]
  }

  # combine metacell objects
  if(length(metacell_list) > 1){
    metacell_obj <- merge(metacell_list[[1]], metacell_list[2:length(metacell_list)])
    
    # need to join layers if this is Seurat 5
    if(CheckSeurat5()){
      metacell_obj <- SeuratObject::JoinLayers(metacell_obj)
    }
  } else{
    metacell_obj <- metacell_list[[1]]
  }

  # set idents for metacell object:
  Idents(metacell_obj) <- metacell_obj@meta.data[[ident.group]]

  # revert to full seurat object if we subsetted earlier
  if(!is.null(cells.use)){
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

  if(CheckSeurat5()){
    param_list['metacell_layer'] <- layer
  } else{
    param_list['metacell_slot'] <- slot
  }

  seurat_obj <- SetWGCNAParams(
    seurat_obj, params = param_list, wgcna_name
  )

  # return the updated seurat object
  seurat_obj
}


#' NormalizeMetacells
#'
#' Wrapper function to run Seurat's NormalizeData function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @export
NormalizeMetacells <- function(seurat_obj, wgcna_name=NULL, ...){
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::NormalizeData(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}

#' ScaleMetacells
#'
#' Wrapper function to run Seurat's ScaleData function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @export
ScaleMetacells <- function(seurat_obj, wgcna_name=NULL, ...){
  if(!exists("features")){
    features = VariableFeatures(seurat_obj)
  }
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::ScaleData(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}

#' RunPCAMetacells
#'
#' Wrapper function to run Seurat's RunPCA function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @export
RunPCAMetacells <- function(seurat_obj, wgcna_name=NULL, ...){
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::RunPCA(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}

#' RunHarmonyMetacells
#'
#' Wrapper function to run harmony's RunHarmony function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @export
RunHarmonyMetacells <- function(seurat_obj, ...){
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- harmony::RunHarmony(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
  seurat_obj
}

#' RunUMAPMetacells
#'
#' Wrapper function to run Seurat's RunUMAP function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
RunUMAPMetacells <- function(seurat_obj, ...){
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::RunUMAP(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
  seurat_obj
}


#' DimPlotMetacells
#'
#' Wrapper function to run Seurat's DimPlot function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
DimPlotMetacells <- function(seurat_obj, ...){
  Seurat::DimPlot(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
}
