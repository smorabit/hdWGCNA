# hdWGCNA 0.1.1.9008 (2022-07-21)
## Added
- New tutorial for consensus co-expression network analysis.

## Changes
- `ModuleNetworkPlot` and `RunModuleUMAP` now checks if `ModuleConnectivity` has been computed in order to throw a more informative error.
- `GetTOM` checks if the TOM file exists in order to throw a more informative error.


# hdWGCNA 0.1.1.9007 (2022-07-20)
## Added
- None

## Changes
- New warning message in `MetacellsByGroups` if there are any groups that are excluded by `min_cells`.
- Assay in Metacell seurat object is now the same as the assay supplied to `MetacellsByGroups`, instead of the default "RNA".
- `ModuleEigengenes` takes `assay` as an argument, clears up some issues with `RunHarmony`.
- `ModuleEigengenes` doesn't require a "counts" slot to be present in the given assay, but now it throws an error if the normalized data slot is missing.

# hdWGCNA 0.1.1.9006 (2022-07-14)
## Added
- None

## Changes
- `ConstructNetwork` now checks if the TOM file already exists, and whether the user wants to overwrite the existing TOM.

# hdWGCNA 0.1.1.9005 (2022-06-17)
## Added
- None

## Changes
- Added new arguments to `MetacellsByGroups` and `ConstructMetacells` to exclude very small groups (`min_cells`), to reach a target number of metacells (`target_metacells`), and to exclude metacells with too much overlap (`max_shared`).

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
