# high dimensional WGCNA

[![R](https://img.shields.io/github/r-package/v/smorabit/hdWGCNA)](https://github.com/smorabit/hdWGCNA/tree/dev)
[![ISSUES](https://img.shields.io/github/issues/smorabit/hdWGCNA)](https://github.com/smorabit/hdWGCNA/issues)
[![Publication](https://img.shields.io/badge/publication-Cell%20Rep%20Meth-%2300A1D7)](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00127-3)
[![Lifecycle:Maturing](https://img.shields.io/badge/Lifecycle-Maturing-007EC6)](https://github.com/smorabit/hdWGCNA)
[![Stars](https://img.shields.io/github/stars/smorabit/hdWGCNA?style=social)](https://github.com/smorabit/hdWGCNA/)

hdWGCNA is an R package for performing weighted gene co-expression
network analysis [(WGCNA)](https://doi.org/10.1186/1471-2105-9-559) in
high dimensional transcriptomics data such as single-cell RNA-seq or
spatial transcriptomics. hdWGCNA is highly modular and can construct
context-specific co-expression networks across cellular and spatial
hierarchies. hdWGNCA identifies modules of highly co-expressed genes and
provides context for these modules via statistical testing and
biological knowledge sources. hdWGCNA uses datasets formatted as
[Seurat](https://satijalab.org/seurat/index.html) objects. Check out the
[hdWGCNA in single-cell data
tutorial](https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html)
or the [hdWGCNA in spatial transcriptomics data
tutorial](https://smorabit.github.io/hdWGCNA/articles/ST_basics.html) to
get started.

**New functionality:** hdWGCNA is now able to perform [Transcription
Factor Regulatory Network
Analysis](https://smorabit.github.io/hdWGCNA/articles/tf_network.html).
This functionality was introduced in our publication [Childs & Morabito
et al., Cell Reports
(2024)](https://www.sciencedirect.com/science/article/pii/S2211124724002845).

If you use hdWGCNA in your research, please cite the manuscript
describing hdWGCNA in addition to the [original WGCNA
publication](https://doi.org/10.1186/1471-2105-9-559):

- [hdWGCNA identifies co-expression networks in high-dimensional
  transcriptomics data. Morabito et al., Cell Reports Methods
  (2023)](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00127-3)

## Installation

We recommend creating an R [conda
environment](https://docs.conda.io/en/latest/) environment for hdWGCNA.

``` bash
# create new conda environment for R
conda create -n hdWGCNA -c conda-forge -c bioconda r-base=4.4 mamba

# activate conda environment
conda activate hdWGCNA

# install critical R packages
mamba install -c conda-forge -c bioconda r-seurat r-hdf5r r-wgcna r-igraph r-tidyverse r-ggraph r-harmony r-enrichr r-devtools

# install Bioconductor packages
mamba install -c conda-forge -c bioconda bioconductor-ucell bioconductor-genomicranges bioconductor-geneoverlap 
```

Next open R and install hdWGCNA.

``` r

# install Bioconductor
install.packages("BiocManager")
BiocManager::install()

# install hdWGCNA from GitHub
devtools::install_github('smorabit/hdWGCNA', ref='dev')
```

Troubleshooting the installation

If you have trouble installing some of the dependencies, we advise
trying to install them using alternative sources. For example, if you
try to install the package using Bioconductor or CRAN, try instead with
conda. Let’s say you have trouble installing `devtools`, try instead to
use conda:

    conda install conda-forge::r-devtools

For specific problems, please write a [GitHub
issue](https://github.com/smorabit/hdWGCNA/issues) and we will try to
help you.

## Suggested Reading

To learn more about the method and to see several examples of analyzing
different datasets, please check out the paper describing hdWGCNA
(Morabito et al 2023). We also suggest our follow-up paper introducing
transcription factor regulatory network analysis with hdWGCNA (Childs &
Morabito et al 2024).

- hdWGCNA method description: [hdWGCNA identifies co-expression networks
  in high-dimensional transcriptomics
  data](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00127-3)
- Transcription factor networks: [Relapse to cocaine seeking is
  regulated by medial habenula NR4A2/NURR1 in
  mice](https://www.sciencedirect.com/science/article/pii/S2211124724002845)

For additional reading, we suggest the original WGCNA publication and
papers describing relevant algorithms for co-expression network
analysis.

- [WGCNA: an R package for weighted correlation network
  analysis](https://doi.org/10.1186/1471-2105-9-559)
- [Defining clusters from a hierarchical cluster tree: the Dynamic Tree
  Cut package for R](https://doi.org/10.1093/bioinformatics/btm563)
- [Is My Network Module Preserved and
  Reproducible?](https://doi.org/10.1371/journal.pcbi.1001057)

**Note about package development:** hdWGCNA is under active development,
so you may run into errors and small typos. We welcome users to write
[GitHub
issues](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue)
to report bugs, ask for help, and to request potential enhancements.
