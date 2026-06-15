# SetupForWGCNA

Create a slot in a Seurat object to store hdWGCNA data

## Usage

``` r
SetupForWGCNA(
  seurat_obj,
  wgcna_name,
  features = NULL,
  metacell_location = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- wgcna_name:

  name of the WGCNA experiment

- features:

  list of features to use for WGCNA

- metacell_location:

  name of the WGCNA experiment to copy the metacell object from

- ...:

  additional parameters to pass to SelectNetworkGenes

## Details

SetupForWGCNA creates a new slot in the Seurat object (seurat_obj) to
store an hdWGCNA experiment with a given name (wgcna_name). This
function calls on SelectNetworkGenes to specify which features will be
used for network analysis. If there is another hdWGCNA experiment
already in the seurat_obj, the same metacell/metaspot object can be used
by specifying the name of the hdWGCNA experiment to the
metacell_location parameter.
