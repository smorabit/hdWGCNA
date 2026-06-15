# PlotSoftPowers

Plot Soft Power Threshold results

## Usage

``` r
PlotSoftPowers(
  seurat_obj,
  selected_power = NULL,
  point_size = 5,
  text_size = 3,
  plot_connectivity = TRUE,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- selected_power:

  power to highlight in the plots

- point_size:

  the size of the points in the plot

- text_size:

  the size of the text in the plot

- plot_connectivity:

  logical indicating whether to plot the connectivity in addition to the
  scale free topplogy fit.

- wgcna_name:

  The name of the WGCNA experiment in seurat_obj

## Examples

``` r
PlotSoftPowers
#> function (seurat_obj, selected_power = NULL, point_size = 5, 
#>     text_size = 3, plot_connectivity = TRUE, wgcna_name = NULL) 
#> {
#>     if (is.null(wgcna_name)) {
#>         wgcna_name <- seurat_obj@misc$active_wgcna
#>     }
#>     pt <- GetPowerTable(seurat_obj, wgcna_name)
#>     if ("group" %in% colnames(pt)) {
#>         power_tables <- pt %>% dplyr::group_split(group)
#>         soft_powers <- sapply(power_tables, function(power_table) {
#>             power_table %>% subset(SFT.R.sq >= 0.8) %>% .$Power %>% 
#>                 min
#>         })
#>     }
#>     else {
#>         if (is.null(selected_power)) {
#>             soft_power <- pt %>% subset(SFT.R.sq >= 0.8) %>% 
#>                 .$Power %>% min
#>         }
#>         else {
#>             soft_power <- selected_power
#>         }
#>         soft_powers <- NULL
#>         power_tables <- list(power = pt)
#>     }
#>     plot_list <- list()
#>     for (i in 1:length(power_tables)) {
#>         pt <- power_tables[[i]]
#>         if (!is.null(soft_powers)) {
#>             soft_power <- soft_powers[i]
#>             print(i)
#>             print(soft_power)
#>             pt <- pt %>% dplyr::select(-group)
#>         }
#>         print(head(pt))
#>         sft_r <- as.numeric(pt[pt$Power == soft_power, "SFT.R.sq"])
#>         mean_k <- as.numeric(pt[pt$Power == soft_power, "mean.k."])
#>         median_k <- as.numeric(pt[pt$Power == soft_power, "median.k."])
#>         max_k <- as.numeric(pt[pt$Power == soft_power, "max.k."])
#>         pt$text_color <- ifelse(pt$Power == soft_power, "white", 
#>             "black")
#>         p1 <- pt %>% ggplot(aes(x = Power, y = SFT.R.sq)) + geom_rect(data = pt[1, 
#>             ], aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.8), 
#>             fill = "grey80", alpha = 0.8, color = NA) + geom_hline(yintercept = sft_r, 
#>             linetype = "dashed") + geom_vline(xintercept = soft_power, 
#>             linetype = "dashed") + geom_point(data = pt[pt$Power == 
#>             soft_power, c("Power", "SFT.R.sq")], aes(x = Power, 
#>             y = SFT.R.sq), inherit.aes = FALSE, color = "black", 
#>             size = point_size) + geom_text(label = pt$Power, 
#>             color = pt$text_color, size = text_size) + scale_y_continuous(limits = c(0, 
#>             1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + ylab("Scale-free Topology Model Fit") + 
#>             xlab("Soft Power Threshold") + theme(axis.line.x = element_blank(), 
#>             axis.line.y = element_blank(), panel.border = element_rect(colour = "black", 
#>                 fill = NA, size = 1))
#>         if (plot_connectivity) {
#>             p2 <- pt %>% ggplot(aes(x = Power, y = mean.k.)) + 
#>                 geom_hline(yintercept = mean_k, linetype = "dashed") + 
#>                 geom_vline(xintercept = soft_power, linetype = "dashed") + 
#>                 geom_point(data = pt[pt$Power == soft_power, 
#>                   c("Power", "mean.k.")], aes(x = Power, y = mean.k.), 
#>                   inherit.aes = FALSE, color = "black", size = point_size) + 
#>                 geom_text(label = pt$Power, color = pt$text_color, 
#>                   size = text_size) + scale_y_continuous(labels = scales::comma) + 
#>                 ylab("Mean Connectivity") + xlab("Soft Power Threshold") + 
#>                 theme(axis.line.x = element_blank(), axis.line.y = element_blank(), 
#>                   panel.border = element_rect(colour = "black", 
#>                     fill = NA, size = 1))
#>             p3 <- pt %>% ggplot(aes(x = Power, y = median.k.)) + 
#>                 geom_hline(yintercept = median_k, linetype = "dashed") + 
#>                 geom_vline(xintercept = soft_power, linetype = "dashed") + 
#>                 geom_point(data = pt[pt$Power == soft_power, 
#>                   c("Power", "median.k.")], aes(x = Power, y = median.k.), 
#>                   inherit.aes = FALSE, color = "black", size = point_size) + 
#>                 geom_text(label = pt$Power, color = pt$text_color, 
#>                   size = text_size) + scale_y_continuous(labels = scales::comma) + 
#>                 ylab("Median Connectivity") + xlab("Soft Power Threshold") + 
#>                 theme(axis.line.x = element_blank(), axis.line.y = element_blank(), 
#>                   panel.border = element_rect(colour = "black", 
#>                     fill = NA, size = 1))
#>             p4 <- pt %>% ggplot(aes(x = Power, y = max.k.)) + 
#>                 geom_hline(yintercept = max_k, linetype = "dashed") + 
#>                 geom_vline(xintercept = soft_power, linetype = "dashed") + 
#>                 geom_point(data = pt[pt$Power == soft_power, 
#>                   c("Power", "max.k.")], aes(x = Power, y = max.k.), 
#>                   inherit.aes = FALSE, color = "black", size = point_size) + 
#>                 geom_text(label = pt$Power, color = pt$text_color, 
#>                   size = text_size) + scale_y_continuous(labels = scales::comma) + 
#>                 ylab("Max Connectivity") + xlab("Soft Power Threshold") + 
#>                 theme(axis.line.x = element_blank(), axis.line.y = element_blank(), 
#>                   panel.border = element_rect(colour = "black", 
#>                     fill = NA, size = 1))
#>             plot_list[[i]] <- list(p1, p2, p3, p4)
#>         }
#>         else {
#>             plot_list[[i]] <- p1
#>         }
#>     }
#>     if (length(plot_list) == 1) {
#>         return(plot_list[[1]])
#>     }
#>     plot_list
#> }
#> <bytecode: 0x55d89ca96578>
#> <environment: namespace:hdWGCNA>
```
