
# single-cell co-expression network analysis <img src="man/figures/logo.png" align="right" height="20%" width="20%" />

[![R](https://img.shields.io/github/r-package/v/smorabit/scWGCNA)](https://github.com/smorabit/scWGCNA/tree/dev)
[![LICENSE](https://img.shields.io/github/license/smorabit/scWGCNA)](LICENSE.md)
[![ISSUES](https://img.shields.io/github/issues/smorabit/scWGCNA)](https://github.com/smorabit/scWGCNA/issues)


scWGCNA is an R package for performing weighted gene co-expression network analysis [(WGCNA)](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) in single-cell
RNA-seq data. scWGCNA constructs co-expression networks in a cell-type-specific manner,
identifies robust modules of inerconnected genes, and provides biological
context for these modules. scWGCNA is directly compatible with
[Seurat](https://satijalab.org/seurat/index.html) objects, one of the most ubiquitous
formats for single-cell data. Check out the [scWGCNA basics tutorial](articles/basic_tutorial.html) to get started.


## Installation

We recommend creating an R [conda environment](https://docs.conda.io/en/latest/)
environment for scWGCNA.

```bash
# create new conda environment for R
conda create -n scWGCNA -c conda-forge r-base r-essentials

# activate conda environment
conda activate scWGCNA
```

Next, open up R and install the required dependencies:

* [Bioconductor](https://www.bioconductor.org/), an R-based software ecosystem for bioinformatics and biostatistics.
* [Seurat](https://satijalab.org/seurat/index.html), a general-purpose toolkit for single-cell data science.
* [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/), a package for co-expression network analysis.
* [igraph](https://igraph.org/r/), a package for general network analysis and visualization.
* [devtools](https://devtools.r-lib.org/), a package for package development in R.

```r
# install BiocManager
install.packages("BiocManager")

# install Bioconductor core packages
BiocManager::install()

# install additional packages:
install.packages(c("Seurat", "WGCNA", "igraph", "devtools"))

```

Now you can install the scWGCNA package using `devtools`.

```r
devtools::install_github('smorabit/scWGCNA', branch='dev')
```
