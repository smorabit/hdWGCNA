# ReassignModules

Reassigns features to different modules

## Usage

``` r
ReassignModules(
  seurat_obj,
  harmonized = TRUE,
  features = NULL,
  new_modules = NULL,
  ignore = FALSE,
  auto_reassign = FALSE,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- harmonized:

  logical indicating whether or not to use harmonized MEs

- features:

  character vector containing features for manual module reassignment,

- new_modules:

  character vector containing modules to reassign the genes

- ignore:

  logical indicating whether or not to ignore error message about
  reassigning non-grey features

- auto_reassign:

  logical, automatically re-assign selected features to the module with
  the highest kME?

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

## Value

seurat_obj with the an updated modules table for the selected wgcna
experiment

## Details

ReassignModules reassigns features with negative kMEs in their assigned
module to the module that had the highest kME for that feature.
Alternatively, this function can manually assign features to different
modules, which can be helpful if certain genes of interest are assigned
to the grey module. We generally do not advise reassigning features to
new modules if they are in non-grey modules, but the user can do so at
their own risk by setting ignore=TRUE.
