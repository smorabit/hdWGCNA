# ConstructNetwork

Constructs a co-expression network and groups genes into modules given a
Seurat object that has been prepared for network analysis.

## Usage

``` r
ConstructNetwork(
  seurat_obj,
  soft_power = NULL,
  min_power = 3,
  tom_outdir = "TOM",
  tom_name = NULL,
  consensus = FALSE,
  overwrite_tom = FALSE,
  wgcna_name = NULL,
  blocks = NULL,
  maxBlockSize = 30000,
  randomSeed = 12345,
  corType = "pearson",
  consensusQuantile = 0.3,
  networkType = "signed",
  TOMType = "signed",
  TOMDenom = "min",
  scaleTOMs = TRUE,
  calibrationQuantile = 0.95,
  sampleForCalibration = TRUE,
  sampleForCalibrationFactor = 1000,
  useDiskCache = TRUE,
  chunkSize = NULL,
  deepSplit = 4,
  pamStage = FALSE,
  detectCutHeight = 0.995,
  minModuleSize = 50,
  mergeCutHeight = 0.2,
  saveConsensusTOMs = TRUE,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- soft_power:

  the soft power used for network construction. Automatically selected
  by default.

- min_power:

  the smallest soft power to be selected if soft_power=NULL

- tom_outdir:

  path to the directory where the TOM will be written

- tom_name:

  prefix name given to the TOM output file

- consensus:

  flag indicating whether or not to perform Consensus network analysis

- wgcna_name:

  name of the WGCNA experiment

- ...:

  additional parameters passed to blockwiseConsensusModules

## Value

seurat_obj with the co-expression network and gene modules computed for
the selected wgcna experiment

## Details

ConstructNetwork builds a co-expression network and identifies clusters
of highly co-expressed genes (modules) from the metacell or metaspot
expression matrix stored in the Seurat object. Before running this
function, the following functions must be run on the input Seurat
object:

1.  SetupForWGCNA

2.  MetacellsByGroups or MetaspotsByGroups, and NormalizeMetacells

3.  SetDatExpr or SetMultiExpr

4.  TestSoftPowers or TestSoftPowersConsensus

This function can also be used to perform consensus network analysis if
consensus=TRUE is selected. ConstructNetwork calls the WGCNA function
blockwiseConsensusModules to compute the adjacency matrix, topological
overlap matrix, and to run the Dynamic Tree Cut algorithm to identify
gene modules. blockwiseConsensusModules has numerous parameters but here
we have selected default parameters that we have found to provide
reasonable results on a variety of single-cell and spatial
transcriptomic datasets.
