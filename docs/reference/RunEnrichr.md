# RunEnrichr

Run Enrichr gene set enrichment tests on hdWGCNA modules

## Usage

``` r
RunEnrichr(
  seurat_obj,
  dbs = c("GO_Biological_Process_2021", "GO_Cellular_Component_2021",
    "GO_Molecular_Function_2021"),
  max_genes = 100,
  wait = TRUE,
  wait_time = 5,
  wgcna_name = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- dbs:

  character vector of EnrichR databases

- max_genes:

  Max number of genes to include per module, ranked by kME.

- wait:

  logical indicating whether or not to wait some time between sending
  requests to the EnrichR server.

- wait_time:

  the number of seconds to wait between sending requests to the EnrichR
  server. Value must be less than 60.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot
