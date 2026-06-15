# This function prepares the expression data for Consensus WGCNA analysis. It populates the `multiExpr` slot in the Seurat object, which contains a list of expression matrices (one for each consensus group, e.g., dataset, sample, condition).

The function operates in three modes:

1.  **Internal Seurat/Metacell Mode (Default):** Extracts expression
    data directly from the Seurat object (either single-cell or metacell
    data) based on `group_name`.

2.  **Pseudobulk Mode (SummarizedExperiment):** Extracts expression data
    from a provided SummarizedExperiment object (passed to `mat`). This
    is the recommended approach for pseudobulk consensus analysis.

3.  **External Matrix Mode:** Extracts expression data from a provided
    large matrix (passed to `mat`) where row names contain delimited
    group identifiers.

## Usage

``` r
SetMultiExpr(
  seurat_obj,
  group_name,
  use_metacells = TRUE,
  group.by = NULL,
  multi.group.by = NULL,
  multi_groups = NULL,
  assay = NULL,
  slot = "data",
  layer = "data",
  mat = NULL,
  mat_group_delim = 3,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object containing the hdWGCNA experiment.

- group_name:

  A string containing the specific group to analyze (e.g., a specific
  cluster or cell type). This filters the data when using Internal
  Seurat/Metacell Mode.

- use_metacells:

  Logical; if TRUE (default), use the metacell expression matrix stored
  in the hdWGCNA experiment. If FALSE, use the full single-cell
  expression matrix. Ignored if `mat` is provided.

- group.by:

  A string containing the name of a column in the Seurat object with
  cell groups (clusters, cell types, etc). If NULL (default), uses
  Seurat Idents.

- multi.group.by:

  A string containing the name of the column that defines the consensus
  groups (e.g., "dataset", "sample", "condition").

  - If using **Seurat/Metacells**, this must be a column in
    `seurat_obj@meta.data`.

  - If using **Pseudobulk (SE)**, this must be a column in
    `colData(mat)`.

- multi_groups:

  A character vector specifying which groups from `multi.group.by` to
  include. If NULL, all unique groups are used.

- assay:

  The name of the assay in the Seurat object (e.g., "RNA", "SCT").

- slot:

  The name of the slot in the Seurat object (e.g., "counts", "data").
  Used for Seurat v4 compatibility.

- layer:

  The name of the layer in the Seurat object (e.g., "counts", "data") OR
  the name of the assay in the SummarizedExperiment (e.g., "VST",
  "counts") if `mat` is provided.

- mat:

  A Matrix or SummarizedExperiment object containing gene expression
  data.

  - **SummarizedExperiment:** The function extracts the assay specified
    by `layer` and subsets columns based on `multi.group.by` in
    `colData`.

  - **Matrix:** The function assumes row names are delimited (e.g.,
    "Cluster1:SampleA") and splits them using `mat_group_delim`.

- mat_group_delim:

  Character; the delimiter used in the row names of `mat` if `mat` is a
  matrix (default is ":"). Ignored if `mat` is a SummarizedExperiment.

- wgcna_name:

  A string containing the name of the WGCNA slot in `seurat_obj@misc`.
  Default = NULL, which retrieves the currently active WGCNA data.

- ...:

  Additional arguments passed to internal helper functions.

## Value

A Seurat object with the `multiExpr` slot populated in the specified
`wgcna_name` experiment.

## Details

This function automatically aligns the genes in the provided data (`mat`
or `seurat_obj`) with the genes selected for WGCNA (via `SetupForWGCNA`
or `GetWGCNAGenes`). It also runs
[`WGCNA::goodGenesMS`](https://rdrr.io/pkg/WGCNA/man/goodGenesMS.html)
to exclude genes with zero variance or excessive missingness across the
consensus groups.
