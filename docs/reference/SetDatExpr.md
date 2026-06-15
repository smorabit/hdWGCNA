# Set Expression Data (Standard WGCNA)

The function operates in three modes:

1.  **Internal Seurat/Metacell Mode (Default):** Extracts expression
    data directly from the Seurat object (either single-cell or metacell
    data).

2.  **Pseudobulk Mode (SummarizedExperiment):** Extracts expression data
    from a provided SummarizedExperiment object (passed to `mat`). This
    is the recommended approach for pseudobulk analysis.

3.  **External Matrix Mode:** Sets the expression data from a provided
    matrix (passed to `mat`).

## Usage

``` r
SetDatExpr(
  seurat_obj,
  group_name,
  use_metacells = TRUE,
  group.by = NULL,
  multi.group.by = NULL,
  multi_group_name = NULL,
  return_seurat = TRUE,
  assay = NULL,
  slot = "data",
  layer = "data",
  mat = NULL,
  features = NULL,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object containing the hdWGCNA experiment.

- group_name:

  A string containing the group to subset the data by (e.g., a specific
  cluster or cell type). Only used if pulling data from `seurat_obj`.

- use_metacells:

  Logical; if TRUE (default), use the metacell expression matrix. If
  FALSE, use the full single-cell expression matrix. Ignored if `mat` is
  provided.

- group.by:

  A string containing the name of a column in the Seurat object with
  cell groups (clusters, cell types, etc). If NULL (default), uses
  Seurat Idents.

- multi.group.by:

  A string containing the name of a column in the Seurat object with
  groups to subset further (e.g. dataset, sample). Only used if pulling
  data from `seurat_obj`.

- multi_group_name:

  A string or character vector specifying which groups from
  `multi.group.by` to include.

- return_seurat:

  Logical; if TRUE (default), returns the Seurat object with the
  `datExpr` slot populated. If FALSE, returns the data frame of
  expression data.

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
    by `layer` and subsets genes to match `GetWGCNAGenes(seurat_obj)`.

  - **Matrix:** The function assumes columns are genes and rows are
    samples.

- features:

  A character vector of genes to use. If NULL (default), uses the genes
  stored in the hdWGCNA experiment.

- wgcna_name:

  A string containing the name of the WGCNA slot in `seurat_obj@misc`.
  Default = NULL, which retrieves the currently active WGCNA data.

- ...:

  Additional arguments passed to
  [`WGCNA::goodGenes`](https://rdrr.io/pkg/WGCNA/man/goodGenes.html).

## Value

A Seurat object with the `datExpr` slot populated in the specified
`wgcna_name` experiment, or a data frame if `return_seurat = FALSE`.

## Details

This function sets up the expression matrix input for standard
(non-consensus) WGCNA based on the metacell expression matrix, the full
expression matrix, or a provided pseudobulk expression matrix.

This function automatically runs
[`WGCNA::goodGenes`](https://rdrr.io/pkg/WGCNA/man/goodGenes.html) to
exclude genes with zero variance or excessive missingness. If `mat` is a
SummarizedExperiment, the function automatically transposes the assay to
the required (Samples x Genes) format.
