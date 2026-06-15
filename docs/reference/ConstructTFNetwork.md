# ConstructTFNetwork

Construct a network of transcription factors and target genes based on
gene co-expression

## Usage

``` r
ConstructTFNetwork(
  seurat_obj,
  model_params,
  nfold = 5,
  callbacks = NULL,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- model_params:

  a list of model parameters to pass to xgboost

- nfold:

  number of folds for cross validation

- callbacks:

  callback functions for the xgboost model

- wgcna_name:

  name of the WGCNA experiment

## Value

seurat_obj with the TF network added

## Details

ConstructTFNetwork uses motif-gene information to build a directed
network of transcription factors (TFs) and target genes. XGBoost
regression is leveraged to model the expression of each gene based on
its candidate TF regulators. This analysis gives us information about
how important each TF is for predicting each gene, allowing us to
prioritize the most likely regulators of each gene. This process is done
on the gene expression matrix stored with SetDatExpr, which is typically
the hdWGCNA metacell gene expression matrix.
