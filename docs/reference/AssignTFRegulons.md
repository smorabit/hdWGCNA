# AssignTFRegulons

Define the set of likely target genes (Regulons) for each transcrition
factor

## Usage

``` r
AssignTFRegulons(
  seurat_obj,
  strategy = "A",
  reg_thresh = 0.01,
  n_tfs = 10,
  n_genes = 50,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- strategy:

  method for defining regulons, "A", "B", or "C". See Details for more
  info.

- reg_thresh:

  threshold for regulatory score in strategies A, B, and C

- n_tfs:

  for strategy A, the number of top TFs to keep for each gene

- n_genes:

  for strategy B, the number of top target genes to keep for each TF

- wgcna_name:

  name of the WGCNA experiment

## Value

seurat_obj with the TF Regulon information added

## Details

AssignTFRegulons uses the TF network information from ConstructTFNetwork
to define sets of confident TF-gene pairs. A "regulon" is the set of
target genes for a given TF, and this function provides three different
strategies to define TF regulons. Strategy "A" selects the top TFs for
each gene, strategy "B" selects the top genes for each TF, and strategy
"C" retains all TF-gene pairs above a certain regulatory score
(reg_thresh).
