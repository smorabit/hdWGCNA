# hdWGCNA 0.1.1.9004 (2022-6-13)
## Added
- None

## Changes
- Bug fix so `ModuleTraitCorrelation` can be run with a single trait.

# hdWGCNA 0.1.1.9003 (2022-6-11)
## Added
- `GetHubGenes` function to extract the top hub genes from the module assignment table.

## Changes
- `ConstructMetacells` stores run statistics as a table, and has a new option to exclude metacells that overlap with each other.

# hdWGCNA 0.1.1.9002 (2022-5-24)
## Added
- `ModuleEigengenes` checks to make sure the data has been scaled.

## Changes
- None

# hdWGCNA 0.1.1.9001 (2022-05-09)
## Added
- Wrote docstring for `PlotKMEs` so it would actually be included in the package.

## Changes
- None

# hdWGCNA 0.1.1.9000 (2022-05-05)
## Added
- `PlotKMEs` function to visualize the genes in a module ranked by kME.

## Changes
- Changed the name of the package from scWGCNA to hdWGCNA, since we plan to accommodate data types beyond single cell alone.
- Updated tutorial to recommend computing kMEs in a specific cell population.

# scWGCNA 0.0.0.9000 (2021-07-15)

- The original (and now unsupported) version of scWGCNA included released in [Nature Genetics paper](https://doi.org/10.1038/s41588-021-00894-z).
