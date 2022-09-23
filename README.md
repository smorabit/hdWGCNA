
# high dimensional WGCNA <img src="man/figures/logo.png" align="right" height="20%" width="20%" />

[![R](https://img.shields.io/github/r-package/v/smorabit/hdWGCNA)](https://github.com/smorabit/hdWGCNA/tree/dev)
[![ISSUES](https://img.shields.io/github/issues/smorabit/hdWGCNA)](https://github.com/smorabit/hdWGCNA/issues)
[![DOI](https://zenodo.org/badge/286864581.svg)](https://zenodo.org/badge/latestdoi/286864581)
[![Stars](https://img.shields.io/github/stars/smorabit/hdWGCNA?style=social)](https://github.com/smorabit/hdWGCNA/)

hdWGCNA is an R package for performing weighted gene co-expression network analysis [(WGCNA)](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) in high dimensional
data such as single-cell RNA-seq or spatial transcriptomics.
hdWGCNA is highly modular and can construct co-expression networks to facilitate multi-scale analysis
of cellular and spatial hierarchies. hdWGNCA identifies robust modules of inerconnected genes, and
provides biologicalcontext for these modules through various biological knowledge sources.
hdWGCNA requires data formatted as [Seurat](https://satijalab.org/seurat/index.html) objects,
one of the most ubiquitous formats for single-cell data. Check out the [hdWGCNA in single-cell data tutorial](https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html) or the [hdWGCNA in spatial transcriptomics data tutorial](https://smorabit.github.io/hdWGCNA/articles/ST_basics.html) to get started.

**Note:** hdWGCNA is under active development, so you may run into errors and small typos. We welcome users to
write [GitHub issues](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue)
to, report bugs, ask for help and ask for potential enhancements. GitHub issues are
preferred to emailing the authors directly.

## Installation

We recommend creating an R [conda environment](https://docs.conda.io/en/latest/)
environment for hdWGCNA.

```bash
# create new conda environment for R
conda create -n hdWGCNA -c conda-forge r-base r-essentials

# activate conda environment
conda activate hdWGCNA
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

Now you can install the hdWGCNA package using `devtools`.

```r
devtools::install_github('smorabit/hdWGCNA', ref='dev')
```

## Suggested Reading

If you are unfamiliar with WGCNA, we suggest reading the original WGCNA publication:

* [WGCNA: an R package for weighted correlation network analysis](https://doi.org/10.1186/1471-2105-9-559)

There are a number of additional relevant publications for WGCNA and related algorithms
like Dynamic Tree Cut and Module Preservation analysis:

* [Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R](https://doi.org/10.1093/bioinformatics/btm563)
* [Eigengene networks for studying the relationships between co-expression modules](https://doi.org/10.1186/1752-0509-1-54)
* [Geometric Interpretation of Gene Coexpression Network Analysis](https://doi.org/10.1371/journal.pcbi.1000117)
* [Is My Network Module Preserved and Reproducible?](https://doi.org/10.1371/journal.pcbi.1001057)

Our original description of applying WGCNA to single-nucleus RNA-seq data:

* [Single-nucleus chromatin accessibility and transcriptomic characterization of Alzheimer’s disease](https://doi.org/10.1038/s41588-021-00894-z)
