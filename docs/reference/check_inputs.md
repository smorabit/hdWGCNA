# Check inputs

Check inputs prior to running to_pseudobulk

## Usage

``` r
check_inputs(
  input,
  meta = meta,
  replicate_col = "replicate",
  cell_type_col = "cell_type",
  label_col = "label"
)
```

## Arguments

- input:

  a single-cell matrix to be converted, with features (genes) in rows
  and cells in columns. Alternatively, a `Seurat`, `monocole3`, or or
  `SingleCellExperiment` object can be directly input.

- meta:

  the accompanying meta data whereby the rownames match the column names
  of `input`.

- replicate_col:

  the vector in `meta` containing the replicate information. Defaults to
  `replicate`.

- cell_type_col:

  the vector in `meta` containing the cell type information. Defaults to
  `cell_type`.

- label_col:

  the vector in `meta` containing the experimental label. Defaults to
  `label`.

- min_cells:

  the minimum number of cells in a cell type to retain it. Defaults to
  `3`.

- min_reps:

  the minimum number of replicates in a cell type to retain it. Defaults
  to `2`.

- min_features:

  the minimum number of expressing cells (or replicates) for a gene to
  retain it. Defaults to `0`.

## Value

a cleaned up expression matrix and meta data object
