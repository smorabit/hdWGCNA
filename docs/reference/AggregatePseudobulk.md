# AggregatePseudobulk

This function aggregates a gene-by-cell count matrix into
gene-by-pseudobulk counts. Pseudobulk groups are defined by the
interaction of a replicate identifier and a group identifier (for
example: sample_id × cell_type). The function builds a sparse design
matrix that maps cells to pseudobulks, multiplies the counts matrix by
that mapping to obtain aggregated counts, filters pseudobulks with too
few contributing cells, removes genes with zero variance across the kept
pseudobulks, and returns a SummarizedExperiment containing assay(s) and
per-pseudobulk metadata (including nCells, nUMI and nFeatures).

## Usage

``` r
AggregatePseudobulk(
  X,
  meta,
  replicate_col,
  group_col,
  min_cells = 10,
  assay_name = "counts"
)
```

## Arguments

- X:

  matrix or Matrix Gene-by-cell count matrix. Can be a base R matrix or
  a sparse Matrix (from the Matrix package). Columns must be cell
  identifiers that match the row names of `meta`.

- meta:

  data.frame Per-cell metadata. Row names must correspond to column
  names of `X`. Must contain the columns specified by `replicate_col`
  and `group_col`.

- replicate_col:

  character(1) Name of the column in `meta` indicating the biological
  replicate (for example sample or individual). Used as the first
  component of the interaction that defines pseudobulks.

- group_col:

  character(1) Name of the column in `meta` indicating the grouping
  factor (for example cell type, cluster, condition). Used as the second
  component of the interaction that defines pseudobulks.

- min_cells:

  integer(1), optional Minimum number of cells required for a pseudobulk
  to be retained. Pseudobulks with strictly greater than `min_cells`
  contributing cells are kept. Default value is 10.

- assay_name:

  character(1), optional Name to assign to the assay in the returned
  SummarizedExperiment. Default is "counts".

## Value

SummarizedExperiment An object with:

- assays: a named list with a single matrix-like assay (genes ×
  pseudobulks) containing aggregated counts. The assay name equals
  `assay_name`.

- colData: a data.frame with one row per pseudobulk (metadata built by
  `make_pseudobulk_metadata(meta, pb_groups)` and subset to kept
  pseudobulks). Additional columns added to colData:

  - nCells: number of cells that contributed to each pseudobulk

  - nUMI: sum of counts across genes for each pseudobulk

  - nFeatures: number of genes with non-zero counts in the pseudobulk

## Details

Create pseudobulk samples by aggregating single-cell (or single-nucleus)
counts according to replicate and group annotations, and return a
SummarizedExperiment.

- Input checks:

  - `X` must be matrix-like (dense matrix or Matrix sparse object).

  - `meta` must be a data.frame with rownames matching `colnames(X)`.

  - `replicate_col` and `group_col` must exist in `meta` and contain no
    NAs.

- Grouping:

  - Pseudobulk groups are created by interaction(meta\[replicate_col\],
    meta\[group_col\], drop = TRUE). This produces factor levels
    representing unique replicate × group combinations.

  - A sparse model matrix (~0 + pb_groups) is constructed to map cells
    to pseudobulks. Columns of the resulting aggregated matrix are
    renamed by removing the "pb_groups" prefix that is created by the
    model matrix.

- Filtering:

  - The function computes the number of cells per pseudobulk (`n_cells`)
    and keeps only pseudobulks with n_cells \> min_cells (strict
    inequality).

  - Genes with zero standard deviation across the retained pseudobulks
    are removed.

- Post-processing:

  - The returned SummarizedExperiment's colData receives nCells, nUMI
    and nFeatures. nUMI is computed as column sums of the aggregated
    counts. nFeatures is computed after thresholding counts to binary
    (counts \> 1 set to 1) and summing per column.

## Edge cases and warnings

- If column names of `X` are not present in rownames(meta), the function
  will stop and report the mismatch (it reports the missing cell ids).

- If all pseudobulks are filtered out by the `min_cells` threshold, the
  function will produce an empty SummarizedExperiment or fail in
  downstream steps; callers should check the returned object.

- The function assumes the existence of a helper
  `make_pseudobulk_metadata()` in the calling environment or package;
  this function must accept the same `meta` and `pb_groups` and return
  rownames corresponding to pseudobulk column order.

## See also

make_pseudobulk_metadata, SummarizedExperiment,
Matrix::sparse.model.matrix
