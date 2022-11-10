
#' ConstructMetaspots
#'
#' Computes metaspots in a given Seurat object containing spatial transcriptomics data.
#' This function is called by MetaspotsByGroups and should NOT be run directly!
#' @param cur_seurat A Seurat object
#' @param mode "sum" or "average"
#' @param assay Assay to extract data for aggregation. Default = 'Spatial'
#' @param slot Slot to extract data for aggregation. Default = 'counts'
#' @keywords ST
#' @export
#' @examples
#'
ConstructMetaspots <- function(
  cur_seurat,
  mode = 'sum',
  assay = 'Spatial',
  slot = 'counts'
){
  # get expression matrix:
  X <- GetAssayData(cur_seurat, slot='counts')

  # check to make sure this is just one sample:
  if(sum(unlist(lapply(Images(cur_seurat), function(x){nrow(cur_seurat@images[[x]]@coordinates) != 0}))) != 1){
    stop("More than one sample present in grouping. Please specify a metadata column with group.by indicating different ST samples.")
  }

  # boundaries
  row_range <- range(cur_seurat$row)
  col_range <- range(cur_seurat$col)

  # loop boundaries
  col_bounds <- col_range[1]:col_range[2]

  if(length(col_bounds) >= 4){
    col_bounds <- col_bounds[which(1:length(col_bounds) %% 4 == 0)]
  } else if(length(col_bounds) == 3){
    warning('The selected grouping is spatially constrained and might fail the metaspot aggregation step. You may wish to change the group.by parameter to form larger groups.')
    col_bounds <- col_bounds[2]
  } else if(length(col_bounds) >= 1){
    warning('The selected grouping is spatially constrained and might fail the metaspot aggregation step. You may wish to change the group.by parameter to form larger groups.')
    col_bounds <- col_bounds[1]
  } else{
    warning('No columns specified. Need to change the group.by parameter.')
    return()
  }

  # even or odd cols?
  if(all(col_bounds %% 2 == 0)){even = TRUE} else{even = FALSE}

  row_bounds <- row_range[1]:row_range[2]
  if(even){
    row_bounds <- row_bounds[row_bounds %% 2 == 0]
  } else{
    row_bounds <- row_bounds[row_bounds %% 2 != 0]
  }

  # note: even number rows go to even number cols etc
  coords <- cur_seurat@meta.data[,c('row', 'col')]

  # compute distances
  distances <- proxy::dist(coords, coords, method='euclidean')

  tmp <- distances[distances != 0]
  min1 <- min(tmp)

  bcs <- c()
  unique_cols <- unique(cur_seurat$col);
  unique_cols <- unique_cols[order(unique_cols)]
  combine_list <- list()

  tmp <- lapply(1:length(col_bounds), function(i){
    x = col_bounds[i]
    tmp <- lapply(1:length(row_bounds), function(j){
      y = row_bounds[j]
      cur_coords <- coords %>% subset(row == y & col == x)
      out <- c()
      if(nrow(cur_coords) > 0){

        cur_bc <- rownames(cur_coords)
        bcs <- c(bcs, cur_bc)

        # get cur distances:
        cur_dist <- distances[cur_bc,]
        cur_dist <- cur_dist[cur_dist != 0]
        cur_dist <- cur_dist[cur_dist == min1]

        # find neighbors that are close to this spot
        ix <- which(unique_cols == x)
        other_bcs <- coords %>% subset(col %in% unique_cols[c(ix-2, ix+2)] & row == y) %>% rownames
        cur_neighbors <- c(names(cur_dist), other_bcs)

        # update list of barcodes:
        combine_list[[cur_bc]] <- cur_neighbors
        if(length(cur_neighbors) < 2){
          return(c())
        }

        # aggregate expression profile for these spots:
        cur_X <- X[,c(cur_bc, cur_neighbors)]
        cur_X <- Matrix::rowSums(cur_X)

        if(mode == 'average'){
          cur_X <- cur_X / (length(cur_neighbors) + 1)
        }

        out <- cur_X

      } else{
        return()
      }
      list(cur_bc, cur_neighbors, out)

    })

    if(length(tmp) <= 1){
      next
    }

    # get bc-neighbor results
    bcs <- unlist(lapply(1:length(tmp), function(k){tmp[[k]][[1]]}))
    cur_neighbors <- lapply(1:length(tmp), function(k){tmp[[k]][[2]]})
    cur_neighbors[sapply(cur_neighbors, is.null)] <- NULL
    names(cur_neighbors) <- bcs

    # combine expression results
    cur_X <- do.call(rbind, lapply(1:length(tmp), function(k){tmp[[k]][[3]]}))

    list(bcs, cur_neighbors, cur_X)

  })

  # get bc-neighbor results
  bcs <- unlist(lapply(1:length(tmp), function(k){tmp[[k]][[1]]}))

  cur_neighbors <- list()
  for(k in 1:length(tmp)){
    cur_neighbors <- c(cur_neighbors, tmp[[k]][[2]])
  }
  names(cur_neighbors) <- bcs

  # combine expression results
  if(length(tmp) <= 1){
    warning('Metaspot aggregation failed for this grouping.')
    return(NULL)
  }
  agg_X <- do.call(rbind, lapply(1:length(tmp), function(k){tmp[[k]][[3]]}))

  # cast to sparse matrix
  agg_X <- as(agg_X, 'sparseMatrix')

  # transpose expression matrix:
  agg_X <- Matrix::t(agg_X)
  colnames(agg_X) <- bcs

  # get metadata:
  cur_meta <- cur_seurat@meta.data[bcs,]

  # make metaspot object:
  metaspot_obj <- CreateSeuratObject(
    counts = agg_X,
    meta = cur_meta,
    assay = assay
  )
  if(slot == 'scale.data'){
    metaspot_obj <- SeuratObject::SetAssayData(
      metaspot_obj,
      slot=slot,
      assay=assay,
      new.data=as.matrix(agg_X)
    )
  }

  # add neighbors:
  metaspot_obj@misc$spot_neighbors <- cur_neighbors

  # compute sparsity:
  agg_X[agg_X > 0] <- 1
  X[X > 0] <- 1
  density_agg <- sum(Matrix::colSums(agg_X) / (nrow(agg_X)*ncol(agg_X)))
  density_orig <- sum(Matrix::colSums(X) / (nrow(X)*ncol(X)))

  metaspot_obj@misc$density_agg <- density_agg
  metaspot_obj@misc$density_orig <- density_orig

  metaspot_obj

}
