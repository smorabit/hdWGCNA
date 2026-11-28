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

test_that('ModuleCorrelogram draws without error', {
  skip_if_not_installed('corrplot')
  local_null_device()
  expect_no_error(ModuleCorrelogram(sim_wgcna('mes')))
})

test_that('ModuleFeaturePlot returns one plot per module', {
  so <- sim_with_umap()
  plots <- ModuleFeaturePlot(so, features = 'hMEs', order_points = TRUE)
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  expect_length(plots, length(mods))
  expect_true(all(vapply(plots, is_gg, logical(1))))
})

test_that('ModuleCorrNetwork draws without error', {
  skip_if_not_installed('qgraph')
  local_null_device()
  expect_no_error(suppressWarnings(ModuleCorrNetwork(sim_with_umap())))
})

test_that('PlotKMEs and DoHubGeneHeatmap return ggplot objects', {
  so <- sim_wgcna('kmes')
  expect_true(is_gg(PlotKMEs(so, n_hubs = 5)))
  expect_true(is_gg(suppressMessages(DoHubGeneHeatmap(so, n_hubs = 5, n_cells = 50))))
})

test_that('ModuleNetworkPlot writes one PDF per module', {
  so <- sim_wgcna('kmes')
  outdir <- withr::local_tempdir()
  suppressMessages(ModuleNetworkPlot(so, outdir = outdir, n_inner = 5, n_outer = 5, n_conns = 50))
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  expect_setequal(list.files(outdir), paste0(mods, '.pdf'))
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
})
