# SetTFRegulons

SetTFRegulons

## Usage

``` r
SetTFRegulons(seurat_obj, tf_regulons, wgcna_name = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object

- tf_regulons:

  dataframe storing the TF regulon info from AssignTFRegulons

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot
