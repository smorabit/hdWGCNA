
# scWGCNA <img src="man/figures/logo.png" align="right" height="20%" width="20%" />

scWGCNA is an R package for performing [weighted gene co-expression network analysis (WGCNA)](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) in single-cell
RNA-seq data. scWGCNA constructs co-expression networks in a cell-type-specific manner,
identifies robust modules of highly inerconnected genes, and provides biological
context for these modules. For ease of use, scWGCNA is directly compatible with
[Seurat](https://satijalab.org/seurat/index.html) objects, one of the most common
formats for single-cell data.


## Installation

We recommend creating an R [conda environment](https://docs.conda.io/en/latest/)
environment for scWGCNA.

```
conda create -n scWGCNA -c conda-forge r-base r-essentials
```

```
install.packages('WGCNA')
install.packages('igraph')
install.packages('devtools')

# install Seurat, check their website for the most up-to-date instructions
install.packages('Seurat')
```

## Installation

Now you can install the scWGCNA package using `devtools`:

```
devtools::install_github('smorabit/scWGCNA')
```
