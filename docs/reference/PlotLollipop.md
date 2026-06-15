# PlotLollipop

An Internal function for PlotDMEsLollipop to Plotting function for the
results of FindDMEs and FindAllDMEs

## Usage

``` r
PlotLollipop(modules, cur_DMEs, pvalue, avg_log2FC = "avg_log2FC")
```

## Arguments

- pvalue:

  p_value or fdr used for the comparison in the DMEs dataframe

- seurat_obj:

  A Seurat object

- DMEs:

  dataframe output from FindDMEs or FindAllDMEs

- group.by:

  the column name of the selected comparison in the DMEs dataframe

- comparison:

  character vector or a list of character vectors containing the
  comparison in the DMEs dataframe

## Value

A ggplot object

## Examples

``` r
plot_list <- PlotLollipop(modules, DMEs, cur_title= c("Group1_vs_Group2"), pvalue, avg_log2FC = 'avg_log2FC')
#> Error in PlotLollipop(modules, DMEs, cur_title = c("Group1_vs_Group2"),     pvalue, avg_log2FC = "avg_log2FC"): could not find function "PlotLollipop"
```
