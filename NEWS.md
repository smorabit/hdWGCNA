# hdWGCNA 0.3.01 (2024-03-07)
## Added
- None

## Changes
- Bugfix in `ModulePreservation`
- Update to module preservation tutorial and project modules tutorial.

# hdWGCNA 0.3.00 (2024-02-27)
## Added
- First version with support for Seurat v5.
- Updated the network visualization tutorial with a tutorial for making custom networks.

## Changes
- Changed `FindDMEs` and `FindAllDMEs` to perform differential testing with module expression scores.


# hdWGCNA 0.2.27 (2024-01-29)
## Added
- New option to `MetacellsByGroups` to specify `dims`.

## Changes
- Fixed `HubGeneNetworkPlot` to allow selecting specific modules.
- `GetHubGenes` now returns genes in order from highest to lowest kME in each module.

# hdWGCNA 0.2.26 (2023-12-05)
## Added
- None 

## Changes
- New options for changing colors in `EnrichrBarPlot`.
- `ConstructNetwork` naming of temporary files updated.

# hdWGCNA 0.2.25 (2023-11-28)
## Added
- None 

## Changes
- Requires Seurat version 4. Will update to support v5 in the future.

# hdWGCNA 0.2.24 (2023-09-28)
## Added
- None 

## Changes
- Update call to Harmony in ModuleEigengenes function

# hdWGCNA 0.2.23 (2023-09-10)
## Added
- None 

## Changes
- Fixed bug in SetDatExpr

# hdWGCNA 0.2.22 (2023-09-08)
## Added
- New tutorial for hdWGCNA with pseudobulk data, including some new functions like `ConstructPseudobulk`. 

## Changes
- Updated `SetDatExpr` and `SetMultiExpr` to use a pseudobulk expression matrix.

# hdWGCNA 0.2.21 (2023-08-31)
## Added
- Additional checks for wgcna_name in several functions.
- New section to the DME tutorial to show how to run it in a loop for multiple clusters.

## Changes
- None.

# hdWGCNA 0.2.20 (2023-08-17)
## Added
- None.

## Changes
- Dependency for tester pacakge.

# hdWGCNA 0.2.19 (2023-06-13)
## Added
- None.

## Changes
- Updated README to include publication, and fixed several igraph function calls.


# hdWGCNA 0.2.18 (2023-04-14)
## Added
- None.

## Changes
- We noticed on rare occasion that EnrichR would give duplicated results for different modules, so we added a new option in `RunEnrichr` to wait in between sending requests to the EnrichR server (default is 5 seconds).

# hdWGCNA 0.2.17 (2023-03-27)
## Added
- None.

## Changes
- New error checks in `SetupForWGCNA` and `SelectNetworkGenes`

# hdWGCNA 0.2.16 (2023-03-20)
## Added
- `PlotDMEsLollipop` function to visualize differential module eigengenes.

## Changes
- New error checks in `MetacellsByGroups`


# hdWGCNA 0.2.15 (2023-03-02)
## Added
- None

## Changes
- Bugfix to allow `ResetModuleNames` and `ResetModuleColors` to work when a grey module is not present.


# hdWGCNA 0.2.14 (2023-02-14)
## Added
- None

## Changes
- Fixed a bug in `SetDatExpr` that would throw an error when group.by=NULL was selected.


# hdWGCNA 0.2.13 (2023-02-13)
## Added
- None

## Changes
- Fixed a bug in `ReassignModules` that cause some modules assignments to be NA for some genes.

# hdWGCNA 0.2.12 (2023-02-04)
## Added
- Module eigengene dynamics with pseudotime tutorial
- `PlotModuleTrajectory` function.

## Changes
- None.

# hdWGCNA 0.2.11 (2023-01-30)
## Added
- MAS-Seq tutorial (still under construction)

## Changes
- Bugfix in `ModuleConnectivity` and `ReassignModules`.

# hdWGCNA 0.2.1 (2023-01-24)
## Added
- `ReassignModules` function.

## Changes
- New option in `ModuleConnectivity` to use `corSparse` to compute the correlation, which greatly reduces runtime and memory usage.
- New option in `ModuleConnectivity` to automatically reassign features to different modules if kME is negative.

# hdWGCNA 0.2.03 (2022-12-15)
## Added
- None

## Changes
- Fixed a bug in `ResetModuleNames`.

# hdWGCNA 0.2.03 (2022-11-10)
## Added
- None

## Changes
- New error checking in `MetaspotsByGroups`
- `MetacellsByGroups` now keeps track of which cells were merged.


# hdWGCNA 0.2.02 (2022-11-01)
## Added
- None

## Changes
- New warning in `SetupForWGCNA` if the user selects a very small number of genes.
- `MetaspotsByGroups` uses sparse matrix format internally.


# hdWGCNA 0.2.01 (2022-10-06)
## Added
- None

## Changes
- Bugfix in `MetaspotsByGroups`.


# hdWGCNA 0.2.00 (2022-09-23)
## Added
- `MetaspotsByGroups` to aggregate neighboring ST spots prior to network analysis.
- Tutorial for applying hdWGCNA to spatial transcriptomics datasets.

## Changes
- None

# hdWGCNA 0.1.2.0001 (2022-09-19)
## Added
- None

## Changes
- networkType option in `TestSoftPowers`.

# hdWGCNA 0.1.2.0000 (2022-09-08)
## Added
- Differential Module Eigengene (DME) tutorial
- FindDMEs function
- FindAllDMEs function
- PlotDMEsVolcano function

## Changes
- None

# hdWGCNA 0.1.1.9015 (2022-09-06)
## Added
- None

## Changes
- New data format check in `MetacellsByGroups` to ensure that the selected slot is present in the selected assay.
- `SetDatExpr` now backs up to the full dataset if the metacell dataset isn't found.
- Changed some text to clarify some points in the Module Preservation tutorial.


# hdWGCNA 0.1.1.9014 (2022-08-26)
## Added
- None

## Changes
- Bugfix in `ModuleConnectivity` that caused kMEs to be out of order.

# hdWGCNA 0.1.1.9013 (2022-08-25)
## Added
- Tutorial for using SCTransform normalized data.

## Changes
- None


# hdWGCNA 0.1.1.9012 (2022-08-24)
## Added
- None

## Changes
- Now includes options for other types of correlations in `ModuleConnectivity`.


# hdWGCNA 0.1.1.9011 (2022-08-17)
## Added
- None

## Changes
- Reverted the `exclude_grey` flag back to doing nothing in the `ModuleEigengenes` function because it messed up some downstream tasks, will resolve in a future update.
- `ProjectModules` now excludes modules with too many missing genes in the query dataset,
tunable by the `overlap_proportion` parameter.


# hdWGCNA 0.1.1.9010 (2022-07-30)
## Added
- None

## Changes
- By default, `ModuleEigengenes` does not compute MEs for the grey module. User can change behavior with the `exclude_grey` flag.


# hdWGCNA 0.1.1.9009 (2022-07-23)
## Added
- None

## Changes
- Bugfix in `ProjectModules`.
- Bugfix in `MetacellsByGroups`

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
