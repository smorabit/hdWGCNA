
#' construct_metacells
#'
#' This function takes a Seurat object and constructs averaged 'metacells' based
#' on neighboring cells.
#' @param seurat_obj A Seurat object
#' @param k Number of nearest neighbors to aggregate. Default = 50
#' @param name A string appended to resulting metalcells. Default = 'agg'
#' @param reduction A dimensionality reduction stored in the Seurat object. Default = 'umap'
#' @param assay Assay to extract data for aggregation. Default = 'RNA'
#' @param slot Slot to extract data for aggregation. Default = 'data'
#' @keywords scRNA-seq
#' @export
#' @examples
#' construct_metacells(pbmc)
construct_metacells <- function(seurat_obj, k=50, name='agg', reduction='umap', assay='RNA', slot='data'){

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

  # make seurat obj:
  seurat_aggr <- CreateSeuratObject(
    counts = new_exprs
  )
  list(seurat_aggr, cell_sample)

}
