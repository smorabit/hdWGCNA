# scWGCNA

scWGCNA is a full-featured bioinformatics package based on [Seurat](https://satijalab.org/seurat/index.html) and [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) to perform co-expression network analysis in single-cell or single-nucleus RNA-seq datasets.
WGCNA was originally built for the analysis of bulk gene expression datasets, and the performance of
vanilla WGCNA on single-cell data is limited due to the inherent sparsity of scRNA-seq data. To account for this,
scWGCNA has a function to aggregate transcriptionally similar cells into pseudo-bulk ***metacells*** before
running the WGCNA pipeline, greatly reducing the sparsity of the dataset while preserving cellular heterogeneity. Furthermore, WGCNA is a well established tool with many different options and parameters,
so we recommend trying different options in network construction that are best suited to your dataset.

## Prerequisites

To run scWGCNA, you first need to have a single-cell transcriptomic dataset in Seurat format with
clustering and dimensionality reduction already computed. If this all sounds like gibberish to you,
I would recommend first looking at the [Seurat guided clustering tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).

## Software requirements

scWGCNA has been tested only on R 3.6 and 4.0 on Mac OS and Linux environments (sorry to all Windows bioinformaticians, if there are any of you out there). To run scWGCNA, there are a few other R packages that you need to install. Open up a R session and enter the following commands:

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

## Citation

If scWGCNA is useful in your research, please consider citing our publication.

* [Single-nucleus chromatin accessibility and transcriptomic characterization of Alzheimer's disease](https://www.nature.com/articles/s41588-021-00894-z)

# scWGCNA tutorial

This section is under construction since I have totally changed how scWGCNA works!
