# ComputeModuleEigengene

Internal helper function that computes module eigengene for a single
module.

## Usage

``` r
ComputeModuleEigengene(
  seurat_obj,
  cur_mod,
  modules,
  group.by.vars = NULL,
  verbose = TRUE,
  vars.to.regress = NULL,
  scale.model.use = "linear",
  pc_dim = 1,
  assay = NULL,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- cur_mod:

  name of a module found in
  seurat_obj@misc\[seurat_obj@misc\$active_wgcna\]\$wgcna_net\$colors

- modules:

  table containing module / gene assignments, as in
  GetModules(seurat_obj).

- group.by.vars:

  groups to harmonize by

- verbose:

  logical indicating whether to print messages

- vars.to.regress:

  character vector of variables in seurat_obj@meta.data to regress when
  running ScaleData

- scale.model.use:

  model to scale data when running ScaleData choices are "linear",
  "poisson", or "negbinom"

- pc_dim:

  Which PC to use as the module eigengene? Default to 1.

- assay:

  Assay in seurat_obj to compute module eigengenes. Default is
  DefaultAssay(seurat_obj)

- wgcna_name:

  name of the WGCNA experiment
