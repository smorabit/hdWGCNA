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
#'   - Pseudobulk groups are created by `interaction(meta[[replicate_col]],
#'     meta[[group_col]], drop = TRUE)`. This produces factor levels representing
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
#' @export
AggregatePseudobulk <- function(
  X,
  meta,
  replicate_col,
  group_col,
  min_cells = 10,
  assay_name = 'counts'
) {
  # ------------------------------------------------------------
  # Input sanity checks
  # ------------------------------------------------------------

  # X must be a matrix-like object
  if (!inherits(X, 'Matrix') && !is.matrix(X)) {
    stop("'X' must be a dense matrix or a sparse 'Matrix' object from the Matrix package.")
  }

  # meta must be a data.frame-like object
  if (!is.data.frame(meta)) {
    stop("'meta' must be a data.frame")
  }

  # Ensure that X and meta have matching cells
  if (!all(colnames(X) %in% rownames(meta))) {
    missing <- setdiff(colnames(X), rownames(meta))
    stop('Mismatch between cells in colnames(X) and rownames(meta)')
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
  pb_groups <- interaction(meta[, replicate_col], meta[, group_col], drop = TRUE)

  # calculate the number of cells per grouping
  n_cells <- table(pb_groups)

  # Create group indicator matrix: cells × pseudobulks
  G <- Matrix::sparse.model.matrix(~ 0 + pb_groups)

  # Pseudobulk = counts × G   (genes × cells  %*%  cells × groups)
  pb <- X %*% G
  colnames(pb) <- gsub('pb_groups', '', colnames(pb))

  # remove pseudobulk reps with insufficient cells
  pb <- pb[, as.logical(n_cells >= min_cells)]

  # calculate the standard deviation of each gene:
  good_genes <- names(which(apply(pb, 1, sd) != 0))
  pb <- pb[good_genes, ]

  # create the pseudobulk meta-data table:
  pb_meta <- make_pseudobulk_metadata(meta, pb_groups)
  pb_meta <- pb_meta[colnames(pb), ]

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
find_replicate_columns <- function(meta, group) {
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
#' @export
NormalizeCounts <- function(
  se,
  method = c('CPM', 'logCPM', 'logNorm', 'VST', 'rlog'),
  assay_name = 'counts',
  new_assay_name = NULL,
  pseudocount = 1,
  ...
) {
  method <- match.arg(method)

  # --- Check inputs ---------------------------------------------------------
  if (!inherits(se, 'SummarizedExperiment')) {
    stop("'se' must be a SummarizedExperiment object.")
  }
  if (!assay_name %in% names(assays(se))) {
    stop(paste0("Assay '", assay_name, "' not found in se."))
  }

  X <- assay(se, assay_name)

  if (!is.matrix(X) && !inherits(X, 'Matrix')) {
    stop('The assay must be a matrix or sparse Matrix.')
  }

  # default name for normalized assay
  if (is.null(new_assay_name)) {
    new_assay_name <- method
  }

  # --- Normalization methods -----------------------------------------------

  ## CPM ---------------------------------------------------------------------
  if (method == 'CPM') {
    lib.size <- colSums(X)
    norm <- t(t(X) / lib.size) * 1e6
  }

  ## log CPM -----------------------------------------------------------------
  if (method == 'logCPM') {
    lib.size <- colSums(X)
    cpm <- t(t(X) / lib.size) * 1e6
    norm <- log2(cpm + pseudocount)
  }

  ## Log-normalization (like Seurat: log1p(counts / size factor * 1e4)) ------
  if (method == 'logNorm') {
    size.factor <- colSums(X) / median(colSums(X))
    scaled <- t(t(X) / size.factor) * 1e4
    norm <- log1p(scaled)
  }

  ## VST using DESeq2 --------------------------------------------------------
  if (method %in% c('VST', 'rlog')) {
    check_suggests('DESeq2', feature = 'NormalizeCounts(method = "VST" or "rlog")')

    # Extract the raw counts and convert to a dense matrix
    mat_dense <- as.matrix(SummarizedExperiment::assay(se, assay_name))

    # Rebuild a SE object with dense assay for DESeq2
    se_dense <- SummarizedExperiment::SummarizedExperiment(
      assays  = list(counts = mat_dense),
      colData = SummarizedExperiment::colData(se)
    )

    # Construct DESeqDataSet
    dds <- DESeq2::DESeqDataSet(se_dense, design = ~1)

    # Run DESeq2 normalization transform
    if (method == 'VST') {
      norm <- SummarizedExperiment::assay(DESeq2::vst(dds, ...))
    } else { # rlog
      norm <- SummarizedExperiment::assay(DESeq2::rlog(dds, ...))
    }
  }

  # --- Add normalized assay and return --------------------------------------
  assay(se, new_assay_name) <- norm
  return(se)
}


#' ConstructPseudobulk (deprecated)
#'
#' This function is deprecated. Use [AggregatePseudobulk()] to build pseudobulk samples,
#' [NormalizeCounts()] to normalize them, and [SetDatExpr()] (with `mat = se`) to use them
#' for network analysis. See the pseudobulk vignette for the full workflow.
#'
#' Constructs a pseudobulk gene expression matrix summarizing the expression of each gene
#' in each biological replicate and group (for example cell type), restricted to the genes
#' selected for the hdWGCNA experiment.
#'
#' @return A matrix with one row per pseudobulk sample (named `group:replicate`) and one
#'   column per hdWGCNA gene.
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing grouping info, ie clusters or celltypes
#' @param replicate_col column in seurat_obj@meta.data denoting each replicate / sample
#' @param label_col Ignored; kept for backwards compatibility.
#' @param assay Assay in seurat_obj to aggregate.
#' @param slot Deprecated, use `layer`.
#' @param layer Layer to extract data for aggregation. Default = 'counts'.
#' @param min_reps The minimum number of different biological replicates allowed. Error will be thrown if the number of reps is too low.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' Pseudobulk samples with fewer than 3 cells are dropped. Genes with the same pseudobulk
#' expression in every sample are dropped by [AggregatePseudobulk()], so the output can
#' contain fewer genes than the hdWGCNA experiment.
#'
#' We recommend at least 20 replicates for pseudobulk network analysis; with fewer
#' replicates, correlations are noisy and results may not be reproducible.
#'
#' @seealso [AggregatePseudobulk()], [NormalizeCounts()], [SetDatExpr()]
#' @keywords internal
#' @export
ConstructPseudobulk <- function(
  seurat_obj,
  group.by,
  replicate_col,
  label_col = NULL,
  assay = 'RNA',
  slot = NULL,
  layer = 'counts',
  min_reps = 20,
  wgcna_name = NULL
) {
  .Deprecated(
    msg = paste0(
      'ConstructPseudobulk() is deprecated. Use AggregatePseudobulk() to build pseudobulk ',
      'samples, NormalizeCounts() to normalize them, and SetDatExpr(mat = se) to use them ',
      'for network analysis.'
    )
  )

  # `slot` is deprecated (defunct in SeuratObject 5); use `layer`
  if (!is.null(slot)) {
    warning('The slot argument is deprecated, please use layer instead.', call. = FALSE)
    layer <- slot
  }

  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # check that selected assay is in the seurat object
  if (!(assay %in% Assays(seurat_obj))) {
    stop(paste0('Invalid choice of assay: ', assay, ' not found in Assays(seurat_obj).'))
  }

  # check that layer is valid
  if (!(layer %in% c('counts', 'data', 'scale.data'))) {
    stop('Invalid choice of layer. Valid choices are counts, data, or scale.data.')
  }

  # check that group.by and replicate_col are valid
  for (col in c(group.by, replicate_col)) {
    if (!(col %in% names(seurat_obj@meta.data))) {
      stop(paste0(col, ' not found in seurat_obj@meta.data.'))
    }
    if (!(is.character(seurat_obj@meta.data[[col]]) || is.factor(seurat_obj@meta.data[[col]]))) {
      stop(paste0(col, ' must be a character or a factor column in seurat_obj@meta.data.'))
    }
  }

  # check how many replicates we have:
  n_reps <- length(unique(seurat_obj@meta.data[[replicate_col]]))
  if (n_reps < min_reps) {
    stop(paste0('The number of biological replicates (', n_reps, ') is smaller than min_reps (', min_reps, ')'))
  }
  if (n_reps < 20) {
    warning(paste0('We strongly recommend at least 20 replicates for pseudobulk hdWGCNA, and there are only ', n_reps, ' replicates detected. Results may not be reproducible or informative with low numbers of replicates so proceed at your own risk.'))
  }

  # aggregate into one gene x replicate matrix per group
  X <- SeuratObject::LayerData(seurat_obj, assay = assay, layer = layer)
  pb_list <- pseudobulk_by_group(
    X,
    meta = seurat_obj@meta.data,
    replicate_col = replicate_col,
    group_col = group.by,
    min_cells = 3
  )

  # merge into one matrix with group:replicate sample names
  datExpr <- do.call(cbind, lapply(names(pb_list), function(x) {
    cur <- pb_list[[x]]
    colnames(cur) <- paste0(x, ':', colnames(cur))
    cur
  }))

  # keep the hdWGCNA genes
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  missing_genes <- setdiff(genes_use, rownames(datExpr))
  if (length(missing_genes) > 0) {
    warning(paste0(
      length(missing_genes), ' hdWGCNA genes have the same pseudobulk expression in every ',
      'sample and were removed.'
    ))
  }
  genes_use <- intersect(genes_use, rownames(datExpr))

  # samples x genes
  t(datExpr[genes_use, , drop = FALSE])
}


# Pseudobulk counts per group: a list with one gene x replicate matrix per level of
# group_col (replicates with fewer than min_cells cells are dropped). Built on
# AggregatePseudobulk(); replaces to_pseudobulk(), which was adapted from Libra.
pseudobulk_by_group <- function(X, meta, replicate_col, group_col, min_cells = 3) {
  se <- AggregatePseudobulk(
    X,
    meta,
    replicate_col = replicate_col,
    group_col = group_col,
    min_cells = min_cells
  )
  pb <- SummarizedExperiment::assay(se)
  pb_meta <- methods::as(SummarizedExperiment::colData(se), 'data.frame')
  groups <- as.character(pb_meta[[group_col]])
  lapply(split(seq_len(ncol(pb)), groups), function(ix) {
    m <- as.matrix(pb[, ix, drop = FALSE])
    colnames(m) <- as.character(pb_meta[[replicate_col]][ix])
    m
  })
}
