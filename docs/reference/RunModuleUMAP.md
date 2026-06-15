# RunModuleUMAP

Run UMAP on co-expression matrix using hub genes as features.

## Usage

``` r
RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,
  exclude_grey = TRUE,
  genes_use = NULL,
  wgcna_name = NULL,
  n_neighbors = 25,
  metric = "cosine",
  spread = 1,
  min_dist = 0.4,
  supervised = FALSE,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- n_hubs:

  number of hub genes to use in the UMAP computation

- exclude_grey:

  logical indicating whether to include grey module

- genes_use:

  character vector of genes to use for the UMAP, must already be present
  in GetModules(seurat_obj)

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- ...:

  Additional parameters supplied to uwot::umap
