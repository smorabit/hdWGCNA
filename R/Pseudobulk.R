#' ConstructPseudobulk
#'
#' Constructs a "pseudobulk" gene expression matrix summarizing the expression levels 
#' of each gene across a grouping variable (cell types for example) in each biological 
#' replicate.
#' 
#' @return a matrix containing pseudobulk expression profiles
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing grouping info, ie clusters or celltypes
#' @param replicate_col column in seurat_obj@meta.data denoting each replicate / sample
#' @param label_col column in seurat_obj@meta.data denoting an additional label of interest, for example disease status or biological sex. This is not a required argument and is typically only used for consensus WGCNA
#' @param assay Assay in seurat_obj containing isoform expression information.
#' @param slot Slot to extract data for aggregation. Default = 'counts'
#' @param layer Layer to extract data for aggregation. Default = 'counts'. Layer is used with Seurat v5 instead of slot.
#' @param min_reps The minimum number of different biological replicates allowed. Error will be thrown if the number of reps is too low. 
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' This function constructs pseudobulk gene expression profiles across the provided cell groups 
#' and the provided biological replicates. We note that low numbers of replicates are typical 
#' in single-cell and spatial transcriptomics due to the large monetary cost of running these experiments, 
#' and pseudobulk-ing your data for hdWGCNA is only recommended in the case where you have a sufficient 
#' number of replicates. Here we have set the minimum recommended number to 20. Using fewer than 20 replicates 
#' risks the results not being reproducible or robust, and therefore are not biologically meaningful due to 
#' spurious noisy correlations.
#' 
#' @import Seurat
#' @import Matrix
#' @export
ConstructPseudobulk <- function(
  seurat_obj,
  group.by,
  replicate_col,
  label_col = NULL,
  assay = 'RNA',
  slot = 'counts',
  layer = 'counts',
  min_reps = 20,
  wgcna_name = NULL
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # check that selected assay is in the seurat object 
  if(!(assay %in% Assays(seurat_obj))){
    stop(paste0('Invalid choice of assay: ', assay, ' not found in Assays(seurat_obj).'))
  }

  # check that slot is valid 
  if(!(slot %in% c('counts', 'data', 'scale.data'))){
    stop('Invalid choice of slot. Valid choices are counts, data, or scale.data.')
  }

  # check that group.by is valid 
  if(!(group.by %in% names(seurat_obj@meta.data))){
    stop(paste0(group.by, ' not found in seurat_obj@meta.data.'))
  }
  if(!(class(seurat_obj@meta.data[,group.by]) %in% c('character', 'factor'))){
    stop('Selected group.by must be a character or a factor, but ', group.by, ' is a ', class(seurat_obj@meta.data[,group.by]), '.')
  }

  # check that replicate_col is valid
  if(!(class(seurat_obj@meta.data[,replicate_col]) %in% c('character', 'factor'))){
    stop('Selected replicate_col must be a character or a factor, but ', replicate_col, ' is a ', class(seurat_obj@meta.data[,replicate_col]), '.')
  }

  if(is.null(label_col)){
    label_col <- replicate_col
  }

  # check that label_col is valid
  if(!(class(seurat_obj@meta.data[,label_col]) %in% c('character', 'factor'))){
    stop('Selected label_col must be a character or a factor, but ', label_col, ' is a ', class(seurat_obj@meta.data[,label_col]), '.')
  }

  # check that the current assay is the selected assay:
  if(DefaultAssay(seurat_obj) != assay){
    stop('DefaultAssay(seurat_obj) is not ', assay, ' please switch the default assay to the desired assay before running this function.')
  }

    # check how many replicates we have:
    n_reps <- length(unique(seurat_obj@meta.data[[replicate_col]]))
    if(n_reps < min_reps){
        stop(paste0('The number of biological replicates (', n_reps, ') is smaller than min_reps (', min_reps, ')'))
    } 
    if(n_reps < 20){
        warning(paste0("We strongly recommend at least 20 replicates for pseudobulk hdWGCNA, and there are only ", n_reps, " replicates detected. Results may not be reproducible or informative with low numbers of replicates so proceed at your own risk."))
    }

    # get the WGCNA genes:
    genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)

    # get expression matrix:
    if(CheckSeurat5()){
      X <- SeuratObject::LayerData(seurat_obj, assay=assay, layer=layer)
    } else{
      X <- Seurat::GetAssayData(seurat_obj, assay=assay, slot=slot)
    }

    # get pseudobulk replicates
    pseudobulk_list <- to_pseudobulk(
        X, meta = seurat_obj@meta.data,
        cell_type_col = group.by,
        replicate_col = replicate_col,
        label_col = label_col,
        min_reps=0
    )

    # merge lists into a matrix
    datExpr <- Reduce(cbind, lapply(names(pseudobulk_list), function(x){
        cur <- pseudobulk_list[[x]]
        colnames(cur) <- paste0(x, ':', colnames(cur))
        cur
    }))

    print('here')
    print(class(datExpr))
    print(length(pseudobulk_list))
    print(dim(pseudobulk_list[[1]]))

    datExpr <- t(datExpr[genes_use,])

    # return the pseudobulk matrix:
    datExpr

}



#' Create a pseudobulk matrix
#' 
#' Convert a single-cell expression matrix (i.e., genes by cells)
#' to a pseudobulk matrix by summarizing counts within biological replicates.
#' This function is 
#' 
#' @param input a single-cell matrix to be converted, with features (genes) in rows
#'   and cells in columns. Alternatively, a \code{Seurat}, \code{monocole3}, or 
#'   or \code{SingleCellExperiment} object can be directly input.
#' @param meta the accompanying meta data whereby the rownames match the column
#'   names of \code{input}.
#' @param replicate_col the vector in \code{meta} containing the replicate 
#'   information. Defaults to \code{replicate}.
#' @param cell_type_col the vector in \code{meta} containing the cell type 
#'   information. Defaults to \code{cell_type}.
#' @param label_col the vector in \code{meta} containing the experimental
#'   label. Defaults to \code{label}. 
#' @param min_cells the minimum number of cells in a cell type to retain it.
#'   Defaults to \code{3}.
#' @param min_reps the minimum number of replicates in a cell type to retain it.
#'   Defaults to \code{2}.
#' @param min_features the minimum number of expressing cells (or replicates) 
#'   for a gene to retain it. Defaults to \code{0}.
#' @return a list of pseudobulk matrices, for each cell type.
#'  
#' @importFrom magrittr %<>% extract
#' @importFrom dplyr %>% rename_ count group_by filter pull n_distinct distinct
#'   summarise
#' @importFrom purrr map map_int
#' @importFrom Matrix rowSums colSums
#' @importFrom stats setNames
to_pseudobulk = function(input, 
                         meta = NULL, 
                         replicate_col = 'replicate',
                         cell_type_col = 'cell_type',
                         label_col = 'label',
                         min_cells = 3,
                         min_reps = 2,
                         min_features = 0,
                         external = T) {
  if (external) {
    # first, make sure inputs are correct
    inputs = check_inputs(
      input, 
      meta = meta,
      replicate_col = replicate_col,
      cell_type_col = cell_type_col,
      label_col = label_col)
    expr = inputs$expr
    meta = inputs$meta
  } else {
    expr = input
  }

  # convert to characters
  meta %<>% mutate(replicate = as.character(replicate),
                   cell_type = as.character(cell_type),
                   label = as.character(label))
  
  # keep only cell types with enough cells
  keep = meta %>%
    dplyr::count(cell_type, label) %>%
    group_by(cell_type) %>%
    dplyr::filter(all(n >= min_cells)) %>%
    pull(cell_type) %>%
    unique()
  
  # process data into gene x replicate x cell_type matrices
  pseudobulks = keep %>%
    map( ~ {
      print(.)
      cell_type = .
      meta0 = meta %>% dplyr::filter(cell_type == !!cell_type)
      expr0 = expr %>% magrittr::extract(, meta0$cell_barcode)
      # catch cell types without replicates or conditions
      if (n_distinct(meta0$label) < 2)
        return(NA)
      replicate_counts = distinct(meta0, label, replicate) %>%
        group_by(label) %>%
        summarise(replicates = n_distinct(replicate)) %>%
        pull(replicates)
      if (any(replicate_counts < min_reps))
        return(NA)
      
      # process data into gene X replicate X cell_type matrice
      mm = model.matrix(~ 0 + replicate:label, data = meta0)
      mat_mm = expr0 %*% mm
      keep_genes = rowSums(mat_mm > 0) >= min_features
      mat_mm = mat_mm[keep_genes, ] %>% as.data.frame()
      mat_mm %<>% as.data.frame()
      colnames(mat_mm) = gsub("replicate|label", "", colnames(mat_mm))
      # drop empty columns
      keep_samples = colSums(mat_mm) > 0
      mat_mm %<>% magrittr::extract(, keep_samples)
      return(mat_mm)
    }) %>%
    setNames(keep)
  
  # drop NAs
  pseudobulks %<>% magrittr::extract(!is.na(.))
  
  # also filter out cell types with no retained genes
  min_dim = map(pseudobulks, as.data.frame) %>% map(nrow)
  pseudobulks %<>% magrittr::extract(min_dim > 1)
  
  # also filter out types without replicates
  min_repl = map_int(pseudobulks, ~ {
    # make sure we have a data frame a not a vector
    tmp = as.data.frame(.)
    targets = data.frame(group_sample = colnames(tmp)) %>%
      mutate(group = gsub(".*\\:", "", group_sample))
    if (n_distinct(targets$group) == 1)
      return(as.integer(0))
    min(table(targets$group))
  })
  pseudobulks %<>% magrittr::extract(min_repl >= min_reps)
  return(pseudobulks)
}

#' Check inputs
#'
#' Check inputs prior to running to_pseudobulk
#'
#'
#' @param input a single-cell matrix to be converted, with features (genes) in rows
#'   and cells in columns. Alternatively, a \code{Seurat}, \code{monocole3}, or
#'   or \code{SingleCellExperiment} object can be directly input.
#' @param meta the accompanying meta data whereby the rownames match the column
#'   names of \code{input}.
#' @param replicate_col the vector in \code{meta} containing the replicate
#'   information. Defaults to \code{replicate}.
#' @param cell_type_col the vector in \code{meta} containing the cell type
#'   information. Defaults to \code{cell_type}.
#' @param label_col the vector in \code{meta} containing the experimental
#'   label. Defaults to \code{label}.
#' @param min_cells the minimum number of cells in a cell type to retain it.
#'   Defaults to \code{3}.
#' @param min_reps the minimum number of replicates in a cell type to retain it.
#'   Defaults to \code{2}.
#' @param min_features the minimum number of expressing cells (or replicates) 
#'   for a gene to retain it. Defaults to \code{0}.
#' @return a cleaned up expression matrix and meta data object
#'
#' @importFrom dplyr %>% rename_ n_distinct mutate_at vars
#' @importFrom magrittr %<>%
#' @importFrom tester is_numeric_matrix is_numeric_dataframe
#' @importFrom methods is
#'
check_inputs = function(input,
                        meta = meta,
                        replicate_col = 'replicate',
                        cell_type_col = 'cell_type',
                        label_col = 'label') {

  # extract cell types and label from metadata
  if ("Seurat" %in% class(input)) {
    # confirm Seurat is installed
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("install \"Seurat\" R package for Augur compatibility with ",
           "input Seurat object", call. = FALSE)
    }
    meta = input@meta.data %>%
      droplevels()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = Seurat::GetAssayData(input, slot = 'counts')
  } else if ("cell_data_set" %in% class(input)) {
    # confirm monocle3 is installed
    if (!requireNamespace("monocle3", quietly = TRUE)) {
      stop("install \"monocle3\" R package for Augur compatibility with ",
           "input monocle3 object", call. = FALSE)
    }
    meta = monocle3::pData(input) %>%
      droplevels() %>%
      as.data.frame()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = monocle3::exprs(input)
  } else if ("SingleCellExperiment" %in% class(input)){
    # confirm SingleCellExperiment is installed
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("install \"SingleCellExperiment\" R package for Augur ",
           "compatibility with input SingleCellExperiment object",
           call. = FALSE)
    }
    meta = SummarizedExperiment::colData(input) %>%
      droplevels() %>%
      as.data.frame()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = SummarizedExperiment::assay(input)
  } else {
    # check if input is sparse matrix or numberic matrix/df
    valid_input = is(input, 'sparseMatrix') ||
      is_numeric_matrix(input) ||
      is_numeric_dataframe(input)
    if (!valid_input)
      stop("input must be Seurat, monocle, sparse matrix, numeric matrix, or ",
           "numeric data frame")
    if (is.null(meta))
      stop("input matrix must be accompanied by a metadata table")
    expr = input
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    labels = as.character(meta[[label_col]])
    cell_types = as.character(meta[[cell_type_col]])
  }
  
  # check dimensions are non-zero
  if (length(dim(expr)) != 2 || !all(dim(expr) > 0)) {
    stop("expression matrix has at least one dimension of size zero")
  }

  # check dimensions match
  n_cells1 = nrow(meta)
  n_cells2 = ncol(expr)
  if (n_cells1 != n_cells2) {
    stop("number of cells in metadata (", n_cells1, ") does not match number ",
         "of cells in expression (", n_cells2, ")")
  }

  # check at least two labels
  if (n_distinct(labels) == 1) {
    stop("only one label provided: ", unique(labels))
  }

  # check for missing labels or cell types
  if (any(is.na(labels))) {
    stop("labels contain ", sum(is.na(labels)), "missing values")
  }
  if (any(is.na(cell_types))) {
    stop("cell types contain ", sum(is.na(cell_types)), "missing values")
  }
  if (!is.null(replicate_col) && any(is.na(replicates))) {
    stop("replicates contain ", sum(is.na(replicates)), "missing values")
  }

  # check for missing replicates
  if (!is.null(replicate_col) && is.null(replicates)) {
    stop("metadata does not contain replicate information")
  }

  # remove missing values
  missing = is.na(expr)
  if (any(missing)) {
    stop("matrix contains ", sum(missing), "missing values")
  }
  
  # clean up the meta data
  if (!is.null(replicate_col)) {
    meta %<>% as.data.frame() %>%
      mutate(cell_barcode = rownames(meta),
             replicate = meta[[replicate_col]],
             cell_type = meta[[cell_type_col]],
             label = meta[[label_col]]) %>%
      mutate_at(vars(replicate, cell_type, label), as.factor)
  } else {
    meta %<>% as.data.frame() %>%
      mutate(cell_barcode = rownames(meta),
             cell_type = meta[[cell_type_col]],
             label = meta[[label_col]]) %>%
      mutate_at(vars(cell_type, label), as.factor)
  }

  # make sure meta contains row names and is a data frame
  rownames(meta) = colnames(expr)
  meta = as.data.frame(meta)
  to_return = list(
    expr = expr,
    meta = meta
  )
  return(to_return)
}