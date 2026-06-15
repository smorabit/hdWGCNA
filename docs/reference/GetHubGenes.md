# GetHubGenes

Extract the top N hub genes for a given set of modules. This function
outputs a table with the gene name, the module, and the kME for that
module for the top N hub genes.

## Usage

``` r
GetHubGenes(seurat_obj, n_hubs = 10, mods = NULL, wgcna_name = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object

- n_hubs:

  the number of hub genes to select for each module

- mods:

  list of modules, selects all modules by default

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot
