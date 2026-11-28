is_gg <- function(x) inherits(x, 'ggplot') || inherits(x, 'patchwork')

test_that('PlotDendrogram returns a dendrogram above the module color bar', {
  so <- sim_wgcna('network')
  p <- PlotDendrogram(so, main = 'test')
  expect_s3_class(p, 'patchwork')
  hc <- GetNetworkData(so)$dendrograms[[1]]
  # three segments per merge, and the color bar follows the leaf order
  expect_equal(nrow(p[[1]]$data), 3 * nrow(hc$merge))
  expect_equal(p[[2]]$data$color, as.character(GetModules(so)$color)[hc$order])
  local_null_device()
  expect_no_error(print(p))
})

test_that('PlotDendrogram still draws the WGCNA version with engine = base', {
  local_null_device()
  expect_null(PlotDendrogram(sim_wgcna('network'), engine = 'base'))
})

test_that('dendrogram_segments places merges like plot.hclust', {
  hc <- stats::hclust(stats::dist(c(1, 2, 6, 7, 20)))
  segs <- dendrogram_segments(hc, hang = -1)
  # leaves start at height 0, every merge reaches its own height
  expect_true(all(segs$y >= 0))
  expect_setequal(unique(segs$yend), hc$height)
  # the root joins the midpoints of its two children
  root <- segs[segs$yend == max(hc$height) & segs$y == segs$yend, ]
  expect_equal(nrow(root), 1)
  expect_true(root$x < root$xend)
})

test_that('ModuleCorrelogram returns a ggplot of module correlations', {
  so <- sim_wgcna('mes')
  mods <- setdiff(colnames(GetMEs(so)), 'grey')
  n <- length(mods)

  p <- suppress_deprecated(ModuleCorrelogram(so))
  expect_s3_class(p, 'ggplot')
  expect_equal(nrow(p$data), n * (n + 1) / 2)
  expect_true(all(p$data$r[p$data$row == p$data$col] == 1))
  expect_true(all(p$data$stars %in% c('', '*', '**', '***')))

  full <- suppress_deprecated(ModuleCorrelogram(so, type = 'full', method = 'tile', order = 'hclust'))
  expect_equal(nrow(full$data), n^2)
  expect_setequal(levels(full$data$col), mods)

  # a second set of features measured in the same cells gives a rectangular matrix
  p2 <- suppress_deprecated(ModuleCorrelogram(so, MEs2 = GetMEs(so, harmonized = FALSE)[, mods[1:2]]))
  expect_equal(nrow(p2$data), 2 * n)
  expect_warning(suppress_deprecated(ModuleCorrelogram(so, tl.cex = 1)), 'no longer uses corrplot')
  expect_warning(ModuleCorrelogram(so), 'deprecated')
  local_null_device()
  expect_no_error(print(p))
})

test_that('ModuleFeaturePlot returns one plot per module', {
  so <- sim_with_umap()
  plots <- ModuleFeaturePlot(so, features = 'hMEs', order_points = TRUE)
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  expect_length(plots, length(mods))
  expect_true(all(vapply(plots, is_gg, logical(1))))
})

test_that('PlotKMEs and DoHubGeneHeatmap return ggplot objects', {
  so <- sim_wgcna('kmes')
  expect_true(is_gg(PlotKMEs(so, n_hubs = 5)))
  expect_true(is_gg(suppressMessages(suppress_deprecated(DoHubGeneHeatmap(so, n_hubs = 5, n_cells = 50)))))
  expect_warning(DoHubGeneHeatmap(so, n_hubs = 5, n_cells = 50), 'deprecated')
})

test_that('ModuleTopologyHeatmap and ModuleTopologyBarplot return ggplot objects', {
  so <- sim_wgcna('kmes')
  mod <- setdiff(levels(GetModules(so)$module), 'grey')[1]
  expect_true(is_gg(suppressMessages(ModuleTopologyHeatmap(so, mod = mod, raster = FALSE))))
  expect_true(is_gg(suppressMessages(ModuleTopologyBarplot(so, mod = mod))))
})

test_that('OverlapModulesDEGs, OverlapDotPlot and OverlapBarPlot work together', {
  skip_if_not_installed('GeneOverlap')
  so <- sim_wgcna('kmes')
  Seurat::Idents(so) <- so$cell_type
  degs <- suppressMessages(Seurat::FindAllMarkers(so, only.pos = TRUE, logfc.threshold = 0.25, verbose = FALSE))
  overlap_df <- OverlapModulesDEGs(so, deg_df = degs, fc_cutoff = 0.25)
  expect_true(all(c('module', 'group', 'odds_ratio', 'pval', 'fdr') %in% colnames(overlap_df)))
  expect_true(is_gg(OverlapDotPlot(overlap_df)))
  bars <- OverlapBarPlot(overlap_df)
  expect_true(all(vapply(bars, is_gg, logical(1))))

  # -log(fdr) with its reference line at -log(0.05), for any plot_var
  fdr_bars <- OverlapBarPlot(overlap_df, plot_var = 'fdr', neglog = TRUE)
  layer_hline <- Filter(function(l) inherits(l$geom, 'GeomHline'), fdr_bars[[1]]$layers)[[1]]
  expect_equal(layer_hline$data$yintercept, -log(0.05))
  expect_true(is_gg(OverlapBarPlot(overlap_df, plot_var = 'pval', neglog = TRUE)[[1]]))
  dots <- OverlapDotPlot(overlap_df, plot_var = 'pval', neglog = TRUE)
  expect_equal(sort(dots$data$pval), sort(-log(overlap_df$pval)))
})

test_that('overlap_transform applies one transform', {
  expect_equal(overlap_transform(c(1, exp(2)), 'x', TRUE, FALSE)$values, c(0, 2))
  expect_equal(overlap_transform(exp(-2), 'p', TRUE, TRUE)$values, 2)
  expect_equal(overlap_transform(5, 'x', FALSE, FALSE)$label, 'x')
})

test_that('PlotKMEs draws without ggplot2 warnings', {
  local_null_device()
  expect_no_warning(print(PlotKMEs(sim_wgcna('kmes'), n_hubs = 5)))
})

test_that('DoHubGeneHeatmap works with a character group column and selected modules', {
  so <- sim_wgcna('kmes')
  so$ct_chr <- as.character(so$cell_type)
  mods <- rev(setdiff(levels(GetModules(so)$module), 'grey'))[1:2]
  p <- suppress_deprecated(DoHubGeneHeatmap(
    so, n_hubs = 5, n_cells = 30, group.by = 'ct_chr', module_names = mods, combine = FALSE
  ))
  expect_length(p, 2)
  expect_error(suppress_deprecated(DoHubGeneHeatmap(so, module_names = 'not_a_module')), 'not_a_module')
})

test_that('ModuleRadarPlot draws one radar per module, scaled to the module maximum', {
  so <- sim_wgcna('mes')
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  plots <- ModuleRadarPlot(so, group.by = 'cell_type', combine = FALSE)
  expect_named(plots, mods)
  # the polygon has one vertex per group; the largest reaches the outer ring
  poly <- plots[[1]]$layers[[3]]$data
  expect_equal(nrow(poly), 3)
  expect_equal(max(sqrt(poly$x^2 + poly$y^2)), 1)
  expect_s3_class(ModuleRadarPlot(so, group.by = 'Sample', ncol = 2), 'patchwork')
  expect_warning(ModuleRadarPlot(so, group.by = 'cell_type', grid.max = 2), 'no longer uses ggradar')
  local_null_device()
  expect_no_error(print(plots[[1]]))
})
