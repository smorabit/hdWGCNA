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
#' @keywords scRNA-seq
#' @export
#' @examples
#' ConstructMetacells(pbmc)
ConstructMetacells <- function(seurat_obj, name='agg', k=50, reduction='umap', assay='RNA', slot='counts',  meta=NULL, return_metacell=FALSE){

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
  while (length(good_choices) > 0 & it < 5000) {
      it <- it + 1
      choice <- sample(seq_len(length(good_choices)), size = 1,
          replace = FALSE)
      new_chosen <- c(chosen, good_choices[choice])
      good_choices <- good_choices[good_choices != good_choices[choice]]
      cell_sample <- nn_map[new_chosen, ]
      others <- seq_len(nrow(cell_sample) - 1)
      this_choice <- cell_sample[nrow(cell_sample), ]
      shared <- sapply(others, get_shared, this_choice = this_choice)

      if (max(shared) < 0.9 * k) {
          chosen <- new_chosen
      }
  }

  cell_sample <- nn_map[chosen, ]
  combs <- combn(nrow(cell_sample), 2)
  shared <- apply(combs, 2, function(x) {
      k2 - length(unique(as.vector(cell_sample[x, ])))
  })

  message(paste0("Overlap QC metrics:\nCells per bin: ",
      k, "\nMaximum shared cells bin-bin: ", max(shared),
      "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ",
      median(shared)))
  if (mean(shared)/k > 0.1)
      warning("On average, more than 10% of cells are shared between paired bins.")

  exprs_old <- GetAssayData(seurat_obj, assay=assay, slot=slot)

  mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in%
      cell_sample[x, , drop = FALSE])
  mask <- Matrix::Matrix(mask)
  new_exprs <- (exprs_old %*% mask) / k
  colnames(new_exprs) <- paste0(name, '_', 1:ncol(new_exprs))
  rownames(cell_sample) <- paste0(name, '_', 1:ncol(new_exprs))
  colnames(cell_sample) <- paste0('knn_', 1:ncol(cell_sample))

  # make seurat obj:
  seurat_aggr <- CreateSeuratObject(
    counts = new_exprs
  )

  # add meta-data:
  if(!is.null(meta)){
    meta_names <- names(meta)
    for(x in meta_names){
      seurat_aggr@meta.data[[x]] <- meta[[x]]
    }
  } else(
    warning('meta not found')
  )

  # add seurat metacell object to the main seurat object:
  if(return_metacell){
    out <- seurat_aggr
  } else{
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- seurat_aggr

    # add other info
    if(is.null(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params)){
      seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params <- list(
        'metacell_k' = k,
        'metacell_reduction' = reduction,
        'metacell_slot' = slot,
        'metacell_assay' = assay
      )
    } else{
      seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[["metacell_k"]] <- k
      seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[["metacell_reduction"]] <- reduction
      seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[["metacell_slot"]] <- slot
      seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[["metacell_assay"]] <- assay
    }

    out <- seurat_obj
  }
  out

}

#' MetacellsByGroups
#'
#' This function takes a Seurat object and constructs averaged 'metacells' based
#' on neighboring cells in provided groupings, such as cluster or cell type.
#' @param seurat_obj A Seurat object
#' @param group.by A character vector of Seurat metadata column names representing groups for which metacells will be computed. Default = 'seurat_clusters'
#' @param k Number of nearest neighbors to aggregate. Default = 50
#' @param name A string appended to resulting metalcells. Default = 'agg'
#' @param reduction A dimensionality reduction stored in the Seurat object. Default = 'umap'
#' @param assay Assay to extract data for aggregation. Default = 'RNA'
#' @param slot Slot to extract data for aggregation. Default = 'data'
#' @keywords scRNA-seq
#' @export
#' @examples
#' MetacellsByGroups(pbmc)
MetacellsByGroups <- function(seurat_obj, group.by=c('seurat_clusters'), k=50, reduction='umap', assay='RNA', slot='counts'){

  # setup grouping variables
  if(length(group.by) > 1){
    seurat_meta <- seurat_obj@meta.data[,group.by]
    for(col in colnames(seurat_meta)){
      seurat_meta[[col]] <- as.character(seurat_meta[[col]])
    }
    seurat_obj$metacell_grouping <- apply(seurat_meta, 1, paste, collapse='_')
  } else {
    seurat_obj$metacell_grouping <- as.character(seurat_obj@meta.data[[group.by]])
  }
  groupings <- unique(seurat_obj$metacell_grouping)

  # unique meta-data for each group
  meta_df <- as.data.frame(do.call(rbind, strsplit(groupings, '_')))
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
    MoreArgs = list(k=k, reduction=reduction, assay=assay, slot=slot, return_metacell=TRUE)
  )
  names(metacell_list) <- groupings

  # combine metacell objects
  metacell_obj <- merge(metacell_list[[1]], metacell_list[2:length(metacell_list)])

  # add seurat metacell object to the main seurat object:
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- metacell_obj

  # add other info
  if(is.null(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params)){
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params <- list(
      'metacell_k' = k,
      'metacell_reduction' = reduction,
      'metacell_slot' = slot,
      'metacell_assay' = assay,
      'metacell_groups' = group.by
    )
  } else{
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[["metacell_k"]] <- k
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[["metacell_reduction"]] <- reduction
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[["metacell_slot"]] <- slot
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[["metacell_assay"]] <- assay
  }
  seurat_obj
}
