# ModuleExprScore

Computes a module score for each co-expression module using Seurat
AddModuleScore or UCell.

## Usage

``` r
ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method = "Seurat",
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- n_genes:

  the number of genes to use for each module, ranked by kME. Setting
  n_genes = 'all' uses all of the genes in a module

- method:

  selected method for module scoring, valid choices are "Seurat" or
  "UCell"

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

seurat_obj with module scores computed for the selected wgcna experiment

## Details

ModuleExprScore provides an alternative function to ModuleEigengenes for
summarizing the expression level of each module. The user can choose
between Seurat AddModuleScore or UCell using the method parameter.
