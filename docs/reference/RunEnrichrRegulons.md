# RunEnrichrRegulons

Run Enrichr gene set enrichment tests on hdWGCNA modules

## Usage

``` r
RunEnrichrRegulons(
  seurat_obj,
  dbs = c("GO_Biological_Process_2021", "GO_Cellular_Component_2021",
    "GO_Molecular_Function_2021"),
  depth = 1,
  use_regulons = TRUE,
  min_genes = 5,
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

- depth:

  Include primary TF target genes (depth=1) or primary + secondary
  (depth=2)

- use_regulons:

  Use the regulons (default=TRUE) or the full set of TF target genes?

- wait:

  logical indicating whether or not to wait some time between sending
  requests to the EnrichR server.

- wait_time:

  the number of seconds to wait between sending requests to the EnrichR
  server. Value must be less than 60.

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot

- max_genes:

  Max number of genes to include per module, ranked by kME.
