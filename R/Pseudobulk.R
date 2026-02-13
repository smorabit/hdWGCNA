
#' AggregatePseudobulk
#'
#' Create pseudobulk samples by aggregating single-cell (or single-nucleus) counts
#' according to replicate and group annotations, and return a SummarizedExperiment.
#'
#' @description
#' This function aggregates a gene-by-cell count matrix into gene-by-pseudobulk
#' counts. Pseudobulk groups are defined by the interaction of a replicate
#' identifier and a group identifier (for example: sample_id × cell_type). The
#' function builds a sparse design matrix that maps cells to pseudobulks,
#' multiplies the counts matrix by that mapping to obtain aggregated counts,
#' filters pseudobulks with too few contributing cells, removes genes with zero
#' variance across the kept pseudobulks, and returns a SummarizedExperiment
#' containing assay(s) and per-pseudobulk metadata (including nCells, nUMI and
#' nFeatures).
#'
#' @param X matrix or Matrix
#'   Gene-by-cell count matrix. Can be a base R matrix or a sparse Matrix
#'   (from the Matrix package). Columns must be cell identifiers that match
#'   the row names of `meta`.
#' @param meta data.frame
#'   Per-cell metadata. Row names must correspond to column names of `X`.
#'   Must contain the columns specified by `replicate_col` and `group_col`.
#' @param replicate_col character(1)
#'   Name of the column in `meta` indicating the biological replicate (for
#'   example sample or individual). Used as the first component of the
#'   interaction that defines pseudobulks.
#' @param group_col character(1)
#'   Name of the column in `meta` indicating the grouping factor (for example
#'   cell type, cluster, condition). Used as the second component of the
#'   interaction that defines pseudobulks.
#' @param min_cells integer(1), optional
#'   Minimum number of cells required for a pseudobulk to be retained. Pseudobulks
#'   with strictly greater than `min_cells` contributing cells are kept. Default
#'   value is 10.
#' @param assay_name character(1), optional
#'   Name to assign to the assay in the returned SummarizedExperiment. Default
#'   is "counts".
#'
#' @return SummarizedExperiment
#'   An object with:
#'   - assays: a named list with a single matrix-like assay (genes × pseudobulks)
#'     containing aggregated counts. The assay name equals `assay_name`.
#'   - colData: a data.frame with one row per pseudobulk (metadata built by
#'     `make_pseudobulk_metadata(meta, pb_groups)` and subset to kept pseudobulks).
#'   Additional columns added to colData:
#'     - nCells: number of cells that contributed to each pseudobulk
#'     - nUMI: sum of counts across genes for each pseudobulk
#'     - nFeatures: number of genes with non-zero counts in the pseudobulk
#'
#' @details
#' - Input checks:
#'   - `X` must be matrix-like (dense matrix or Matrix sparse object).
#'   - `meta` must be a data.frame with rownames matching `colnames(X)`.
#'   - `replicate_col` and `group_col` must exist in `meta` and contain no NAs.
#' - Grouping:
#'   - Pseudobulk groups are created by interaction(meta[[replicate_col]],
#'     meta[[group_col]], drop = TRUE). This produces factor levels representing
#'     unique replicate × group combinations.
#'   - A sparse model matrix (~0 + pb_groups) is constructed to map cells to
#'     pseudobulks. Columns of the resulting aggregated matrix are renamed by
#'     removing the "pb_groups" prefix that is created by the model matrix.
#' - Filtering:
#'   - The function computes the number of cells per pseudobulk (`n_cells`) and
#'     keeps only pseudobulks with n_cells > min_cells (strict inequality).
#'   - Genes with zero standard deviation across the retained pseudobulks are
#'     removed.
#' - Post-processing:
#'   - The returned SummarizedExperiment's colData receives nCells, nUMI and
#'     nFeatures. nUMI is computed as column sums of the aggregated counts.
#'     nFeatures is computed after thresholding counts to binary (counts > 1
#'     set to 1) and summing per column.
#'
#' @section Edge cases and warnings:
#' - If column names of `X` are not present in rownames(meta), the function
#'   will stop and report the mismatch (it reports the missing cell ids).
#' - If all pseudobulks are filtered out by the `min_cells` threshold, the
#'   function will produce an empty SummarizedExperiment or fail in downstream
#'   steps; callers should check the returned object.
#' - The function assumes the existence of a helper `make_pseudobulk_metadata()`
#'   in the calling environment or package; this function must accept the same
#'   `meta` and `pb_groups` and return rownames corresponding to pseudobulk
#'   column order.
#' 
#' @seealso
#' make_pseudobulk_metadata, SummarizedExperiment, Matrix::sparse.model.matrix
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay assay<- assays colData colData<-
#' @importFrom Matrix sparse.model.matrix
#' @export
AggregatePseudobulk <- function(
    X,
    meta,
    replicate_col, 
    group_col,
    min_cells = 10,
    assay_name = 'counts'
){

    # ------------------------------------------------------------
    # Input sanity checks
    # ------------------------------------------------------------

    # X must be a matrix-like object
    if (!inherits(X, "Matrix") && !is.matrix(X)) {
        stop("'X' must be a dense matrix or a sparse 'Matrix' object from the Matrix package.")
    }

    # meta must be a data.frame-like object
    if (!is.data.frame(meta)) {
        stop("'meta' must be a data.frame")
    }

    # Ensure that X and meta have matching cells
    if (!all(colnames(X) %in% rownames(meta))) {
        missing <- setdiff(colnames(X), rownames(meta))
        stop("Mismatch between cells in colnames(X) and rownames(meta)")
    }

    # Optionally reorder meta to match X
    meta <- meta[colnames(X), , drop = FALSE]

    cols_to_check <- c(replicate_col, group_col)

    for (col in cols_to_check) {

        # Ensure column exists
        if (!(col %in% colnames(meta))) {
            stop(sprintf("Column '%s' not found in meta.", col))
        }

        # Check for missing values
        if (any(is.na(meta[[col]]))) {
            stop(sprintf("Column '%s' contains NA values.", col))
        }
    }

    # define the pseudobulk grouping based on replicate_col
    pb_groups <- interaction(meta[,replicate_col], meta[,group_col], drop=TRUE)
    
    # calculate the number of cells per grouping
    n_cells <- table(pb_groups)

    # Create group indicator matrix: cells × pseudobulks
    G <- Matrix::sparse.model.matrix(~0 + pb_groups)

    # Pseudobulk = counts × G   (genes × cells  %*%  cells × groups)
    pb <- X %*% G
    colnames(pb) <- gsub('pb_groups', '', colnames(pb))

    # remove pseudobulk reps with insufficient cells
    pb <- pb[,as.logical(n_cells >= min_cells)]

    # calculate the standard deviation of each gene:
    good_genes <- names(which(apply(pb, 1, sd) != 0))
    pb <- pb[good_genes,]

    # create the pseudobulk meta-data table:
    pb_meta <- make_pseudobulk_metadata(meta, pb_groups)
    pb_meta <- pb_meta[colnames(pb),]

    assay_list <- list('tmp' = pb)
    names(assay_list) <- assay_name

    # Create a SummarizedExperiment object
    se <- SummarizedExperiment(
        assays = assay_list,
        colData = pb_meta
    )

    # add the number of cells
    colData(se)$nCells <- as.numeric(n_cells[colnames(se)])

    # nUMI and nFeatures detected:
    colData(se)$nUMI <- colSums(pb)
    pb[pb > 1] <- 1
    colData(se)$nFeatures <- colSums(pb)

    # return the SummarizedExperiment object:
    return(se)
}

#' Find Replicate Columns
#'
#' Internal helper to identify metadata columns that are invariant within 
#' pseudobulk groups.
#'
#' @param meta data.frame of cell metadata.
#' @param group factor defining the pseudobulk groups.
#'
#' @return Character vector of column names.
#' @noRd
#' @keywords internal
find_replicate_columns <- function(meta, group){
  is_replicate_col <- function(col) {
    # For each group, check if all entries in that cluster-sample group are identical
    all(tapply(col, group, function(x) length(unique(x)) == 1))
  }
  
  replicate_cols <- names(meta)[sapply(meta, is_replicate_col)]
  replicate_cols
}

#' Create Pseudobulk Metadata
#'
#' Internal helper function to collapse single-cell metadata into pseudobulk-level
#' metadata. It iterates through the metadata columns and identifies those that
#' are consistent (invariant) within each pseudobulk group (e.g., "Age", "Sex",
#' "Condition"), retaining them in the output.
#'
#' @param meta data.frame
#'   The original cell-level metadata.
#' @param group factor
#'   A factor vector defining the pseudobulk groups (e.g. interaction of sample and cluster).
#'   Must be the same length as the number of rows in `meta`.
#'
#' @return data.frame
#'   A data frame with one row per pseudobulk group and columns corresponding
#'   to the invariant metadata fields.
#' @noRd
#' @keywords internal
make_pseudobulk_metadata <- function(meta, group) {
  
  # identify which columns are consistent within groups
  replicate_cols <- find_replicate_columns(meta, group)
  unique_levels <- levels(group)
  first_indices <- match(unique_levels, group)
  
  # subset the original metadata using these indices
  out <- meta[first_indices, replicate_cols, drop = FALSE]
  rownames(out) <- unique_levels
  return(out)
}

#' NormalizeCounts
#'
#' Perform per-sample normalization of a count assay stored in a
#' SummarizedExperiment and add the normalized matrix as a new assay.
#'
#' @param se A SummarizedExperiment containing a raw counts assay.
#' @param method Character; one of "CPM", "logCPM", "logNorm", "VST", "rlog".
#'   - "CPM": counts per million.
#'   - "logCPM": log2(CPM + pseudocount).
#'   - "logNorm": Seurat-style log1p(counts / size_factor * 1e4) where
#'     size_factor = colSums(counts) / median(colSums(counts)).
#'   - "VST", "rlog": variance-stabilizing transform or rlog via DESeq2.
#' @param assay_name Character scalar; name of the assay in `se` to normalize
#'   (default: "counts").
#' @param new_assay_name Character or NULL; name to assign the normalized assay.
#'   If NULL, defaults to the chosen `method`.
#' @param pseudocount Numeric scalar added to CPM before log2 in "logCPM"
#'   (default: 1).
#' @param ... Additional arguments forwarded to DESeq2::vst or DESeq2::rlog when
#'   `method` is "VST" or "rlog".
#'
#' @return A SummarizedExperiment identical to `se` but with a new assay named
#'   `new_assay_name` containing the normalized matrix.
#'
#' @details The function checks that `se` is a SummarizedExperiment and that
#'   `assay_name` exists and is a (possibly sparse) matrix. For "VST" and
#'   "rlog", DESeq2 must be installed; the function converts the assay to a
#'   dense matrix and constructs a DESeqDataSet with design ~ 1 before
#'   applying the transform. Errors are raised for invalid inputs.
#'
#' @examples
#' # Basic usage (assuming `se` is a SummarizedExperiment with a "counts" assay)
#' # se_norm <- NormalizeSE(se, method = "logCPM")
#'
#' @seealso DESeq2::vst, DESeq2::rlog
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay assay<- assays colData colData<-
#' @importFrom Matrix sparse.model.matrix
#' @export
NormalizeCounts <- function(
    se,
    method = c("CPM", "logCPM", "logNorm", "VST", "rlog"),
    assay_name = "counts",
    new_assay_name = NULL,
    pseudocount = 1,
    ...
){
    method <- match.arg(method)

    # --- Check inputs ---------------------------------------------------------
    if (!inherits(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object.")
    }
    if (!assay_name %in% names(assays(se))) {
        stop(paste0("Assay '", assay_name, "' not found in se."))
    }

    X <- assay(se, assay_name)

    if (!is.matrix(X) && !inherits(X, "Matrix")) {
        stop("The assay must be a matrix or sparse Matrix.")
    }

    # default name for normalized assay
    if (is.null(new_assay_name)) {
        new_assay_name <- method
    }

    # --- Normalization methods -----------------------------------------------
    
    ## CPM ---------------------------------------------------------------------
    if (method == "CPM") {
        lib.size <- colSums(X)
        norm <- t(t(X) / lib.size) * 1e6
    }

    ## log CPM -----------------------------------------------------------------
    if (method == "logCPM") {
        lib.size <- colSums(X)
        cpm <- t(t(X) / lib.size) * 1e6
        norm <- log2(cpm + pseudocount)
    }

    ## Log-normalization (like Seurat: log1p(counts / size factor * 1e4)) ------
    if (method == "logNorm") {
        size.factor <- colSums(X) / median(colSums(X))
        scaled <- t(t(X) / size.factor) * 1e4
        norm <- log1p(scaled)
    }

    ## VST using DESeq2 --------------------------------------------------------
    if (method %in% c("VST", "rlog")) {
        if (!requireNamespace("DESeq2", quietly = TRUE)) {
            stop("DESeq2 must be installed for VST / rlog.")
        }

        # Extract the raw counts and convert to a dense matrix
        mat_dense <- as.matrix(SummarizedExperiment::assay(se, assay_name))

        # Rebuild a SE object with dense assay for DESeq2
        se_dense <- SummarizedExperiment::SummarizedExperiment(
            assays  = list(counts = mat_dense),
            colData = SummarizedExperiment::colData(se)
        )

        # Construct DESeqDataSet
        dds <- DESeq2::DESeqDataSet(se_dense, design = ~ 1)

        # Run DESeq2 normalization transform
        if (method == "VST") {
            norm <- SummarizedExperiment::assay(DESeq2::vst(dds, ...))
        } else {  # rlog
            norm <- SummarizedExperiment::assay(DESeq2::rlog(dds, ...))
        }
    }

    # --- Add normalized assay and return --------------------------------------
    assay(se, new_assay_name) <- norm
    return(se)
}









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
#' @keywords internal
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