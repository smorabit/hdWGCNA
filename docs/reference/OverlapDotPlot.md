# OverlapDotPlot

Makes barplots from Enrichr data

## Usage

``` r
OverlapDotPlot(
  overlap_df,
  plot_var = "odds_ratio",
  logscale = TRUE,
  neglog = FALSE,
  plot_significance = TRUE,
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

- plot_significance:

  logical controlling whether to plot the significance levels on top of
  the dots

## Examples

``` r
OverlapDotPlot
#> function (overlap_df, plot_var = "odds_ratio", logscale = TRUE, 
#>     neglog = FALSE, plot_significance = TRUE, ...) 
#> {
#>     label <- plot_var
#>     if (logscale) {
#>         overlap_df[[plot_var]] <- log(overlap_df[[plot_var]])
#>         label <- paste0("log(", plot_var, ")")
#>     }
#>     if (neglog) {
#>         overlap_df[[plot_var]] <- -1 * overlap_df[[plot_var]]
#>         label <- paste0("-", label)
#>     }
#>     p <- overlap_df %>% ggplot(aes(x = module, y = group)) + 
#>         geom_point(aes(size = get(plot_var), alpha = get(plot_var)), 
#>             color = overlap_df$color) + RotatedAxis() + ylab("") + 
#>         xlab("") + labs(size = label, alpha = label) + theme(plot.title = element_text(hjust = 0.5), 
#>         axis.line.x = element_blank(), axis.line.y = element_blank(), 
#>         panel.border = element_rect(colour = "black", fill = NA, 
#>             size = 1))
#>     if (plot_significance) {
#>         p <- p + geom_text(aes(label = Significance))
#>     }
#>     p
#> }
#> <bytecode: 0x55d89e135788>
#> <environment: namespace:hdWGCNA>
```
