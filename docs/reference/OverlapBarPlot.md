# OverlapBarPlot

Plots the results from OverlapModulesDEGs as a bar plot

## Usage

``` r
OverlapBarPlot(
  overlap_df,
  plot_var = "odds_ratio",
  logscale = FALSE,
  neglog = FALSE,
  label_size = 2,
  ...
)
```

## Arguments

- overlap_df:

  the Module/DEG overlap table from OverlapModulesDEGs

- plot_var:

  the name of the overlap statistic to plot

- logscale:

  logical controlling whether to plot the result on a log scale, useful
  for odds ratio

- neglog:

  logical controlling wehether to plot the result as a negative log,
  useful for p-value / FDR

- label_size:

  the size of the module labels in the bar plot

## Examples

``` r
OverlapBarPlot
#> function (overlap_df, plot_var = "odds_ratio", logscale = FALSE, 
#>     neglog = FALSE, label_size = 2, ...) 
#> {
#>     label <- plot_var
#>     if (plot_var == "odds_ratio") {
#>         yint <- 1
#>     }
#>     else if (plot_var == "fdr") {
#>         yint <- 0.05
#>     }
#>     if (logscale) {
#>         overlap_df[[plot_var]] <- log(overlap_df[[plot_var]])
#>         label <- paste0("log(", plot_var, ")")
#>         yint = log(yint)
#>     }
#>     if (neglog) {
#>         overlap_df[[plot_var]] <- -1 * log(overlap_df[[plot_var]])
#>         label <- paste0("-log(", label, ")")
#>         yint = -1 * log(yint)
#>     }
#>     groups <- overlap_df$group %>% as.character %>% unique
#>     plot_list <- list()
#>     for (cur_group in groups) {
#>         cur_df <- overlap_df %>% subset(group == cur_group)
#>         p <- cur_df %>% ggplot(aes(x = reorder(module, get(plot_var)), 
#>             y = get(plot_var))) + geom_bar(stat = "identity", 
#>             fill = cur_df$color) + coord_flip() + xlab("") + 
#>             ylab(label) + ggtitle(cur_group) + theme(axis.line.y = element_blank(), 
#>             axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
#>             plot.title = element_text(hjust = 0.5))
#>         if (plot_var == "fdr" | plot_var == "odds_ratio") {
#>             p <- p + geom_hline(yintercept = yint, linetype = "dashed", 
#>                 color = "gray")
#>         }
#>         p <- p + geom_text(aes(label = module, x = module, y = get(plot_var)), 
#>             color = "black", size = label_size, hjust = "inward")
#>         plot_list[[cur_group]] <- p
#>     }
#>     plot_list
#> }
#> <bytecode: 0x55d89b0e7560>
#> <environment: namespace:hdWGCNA>
```
