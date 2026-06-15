# SelectNetworkGenes

Select genes that will be used for co-expression network analysis

## Usage

``` r
SelectNetworkGenes(
  seurat_obj,
  gene_select = "variable",
  fraction = 0.05,
  group.by = NULL,
  gene_list = NULL,
  assay = NULL,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- fraction:

  A numeric that determines the minimum cells that a gene must be
  expressed in order to be included. For example, fraction = 0.05 means
  that 5% of cells must express a gene (count \> 0) for it to be
  included.

- gene_list:

  A character string of gene names, only used if type = "custom"

- assay:

  Assay in seurat_obj to compute module eigengenes. Default is
  DefaultAssay(seurat_obj)

- wgcna_name:

  name of the WGCNA experiment

- type:

  How to select genes? Select "variable", "fraction", "all", or
  "custom".

## Details

SelectNetworkGenes allows us to specify the genes that will be used for
co-expression network analysis. This function is called by
SetupForWGCNA. By default, the variable features in
VariableFeatures(seurat_obj) are used. A custom gene list can also be
used with the gene_list parameter and setting gene_select="custom".

We can also identify genes that are expressed above 0 in a certain
proportion of the dataset by settting gene_select='fraction'. For
example, by setting fraction=0.1 and group.by='seurat_clusters', this
function will identify the set of genes that are expressed in 10% of
cells in at least one of the clusters.
