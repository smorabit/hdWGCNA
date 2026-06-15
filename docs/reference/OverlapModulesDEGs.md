# OverlapModulesDEGs

Performs Fisher's Exact Test for overlap between DEGs and hdWGCNA
modules.

## Usage

``` r
OverlapModulesDEGs(
  seurat_obj,
  deg_df,
  fc_cutoff = 0.5,
  group_col = "cluster",
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- deg_df:

  DEG table formatted like the output from Seurat's FindMarkers

- fc_cutoff:

  log fold change cutoff for DEGs to be included in the overlap test

- group_col:

  the name of the column in deg_df containing the cell grouping
  information

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot
