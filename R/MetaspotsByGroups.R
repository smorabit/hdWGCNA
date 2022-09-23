
#' MetaspotsByGroups
#'
#' Computes metaspots in a given Seurat object containing spatial transcriptomics data.
#'
#' @param cur_seurat A Seurat object
#' @param group.by A character vector of Seurat metadata column names representing groups for which metacells will be computed.
#' @param assay Assay to extract data for aggregation. Default = 'Spatial'
#' @param slot Slot to extract data for aggregation. Default = 'counts'
#' @param mode determines how to make gene expression profiles for metacells from their constituent single cells. Options are "average" or "sum".
#' @param wgcna_name name of the WGCNA experiment
#' @keywords ST
#' @export
#' @examples
#'
MetaspotsByGroups <- function(
  seurat_obj,
  group.by=c('seurat_clusters'),
  ident.group='seurat_clusters',
  assay = 'Spatial',
  slot = 'counts',
  mode = 'sum',
  wgcna_name = NULL
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

  # check that the image coordinates are present:
  if(!all(c('row', 'col', 'imagerow', 'imagecol') %in% colnames(seurat_obj@meta.data))){
    stop('Spatial coordinates missing from seurat_obj@meta.data, must have columns named row, col, imagerow, and imagecol.')
  }

  # check assay:
  if(is.null(assay)){
    assay <- DefaultAssay(seurat_obj)
  } else if(!(assay %in% names(seurat_obj@assays))){
    stop(paste0('Assay ', assay, ' not found in seurat_obj. Select a valid assay: ', paste0(names(seurat_obj@assays), collapse = ', ')))
  }

  # check slot:
  if(!(slot %in% c('counts', 'data', 'scale.data'))){
    stop('Invalid input for slot. Valid choices are counts, data, scale.data.')
  } else{

    # check the shape of the slot
    slot_dim <- dim(GetAssayData(seurat_obj, assay=assay, slot=slot))
    if(any(slot_dim) == 0){
      stop(paste(c("Selected slot ", slot, " not found in this assay.")))
    }
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

  # split seurat obj by groupings
  seurat_list <- lapply(groupings, function(x){seurat_obj[,seurat_obj$metacell_grouping == x]})
  names(seurat_list) <- groupings

  # compute ConstructMetaspots
  metaspot_list <- mapply(
    ConstructMetaspots,
    cur_seurat = seurat_list,
    MoreArgs = list(mode=mode, assay=assay, slot=slot)
  )
  names(metaspot_list) <- groupings


  # remove NULL
  remove <- which(sapply(metaspot_list, is.null))
  if(length(remove) > 1){
    metaspot_list <- metaspot_list[-remove]
  }

  # combine density information:
  if(length(metaspot_list) > 1){
    density_agg <- unlist(lapply(metaspot_list, function(x){x@misc$density_agg}))
    density_orig <- unlist(lapply(metaspot_list, function(x){x@misc$density_orig}))
    density_names <- names(density_agg)
  } else{
    density_agg <- metaspot_list[[1]]@misc$density_agg
    density_orig <- metaspot_list[[1]]@misc$density_orig
    density_names <- 'all'
  }

  density_df <- data.frame(
    group = density_names,
    agg = as.numeric(density_agg),
    orig = as.numeric(density_orig)
  )

  spot_neighbors <- list()
  for(k in 1:length(metaspot_list)){
    spot_neighbors <- c(spot_neighbors, metaspot_list[[k]]@misc$spot_neighbors)
  }

  # combine metacell objects
  if(length(metaspot_list) > 1){
    metaspot_obj <- merge(metaspot_list[[1]], metaspot_list[2:length(metaspot_list)])
  } else{
    metaspot_obj <- metaspot_list[[1]]
  }

  # add density and neighbor info
  metaspot_obj@misc$spot_neighbors <- spot_neighbors
  metaspot_obj@misc$density <- density_df

  # set idents for metaspot object:
  Idents(metaspot_obj) <- metaspot_obj@meta.data[[ident.group]]

  # add metaspot object to dataset
  seurat_obj <- SetMetacellObject(seurat_obj, metaspot_obj, wgcna_name)

  seurat_obj

}
