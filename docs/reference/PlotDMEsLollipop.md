# PlotDMEsLollipop

Plotting function for the results of FindDMEs and FindAllDMEs

## Usage

``` r
PlotDMEsLollipop(
  seurat_obj,
  DMEs,
  wgcna_name,
  group.by = NULL,
  comparison = NULL,
  pvalue,
  avg_log2FC = "avg_log2FC"
)
```

## Arguments

- seurat_obj:

  A Seurat object

- DMEs:

  dataframe output from FindDMEs or FindAllDMEs

- group.by:

  the column name of the selected comparison in the DMEs dataframe

- comparison:

  character vector or a list of character vectors containing the
  comparison in the DMEs dataframe

- pvalue:

  p_value or fdr used for the comparison in the DMEs dataframe

## Value

A ggplot object PlotDMEsLollipop
