# SetTFEval

SetTFEval

## Usage

``` r
SetTFEval(seurat_obj, tf_eval, wgcna_name = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object

- tf_eval:

  dataframe storing the TF network evaluation info from
  ConstructTFNetwork

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot
