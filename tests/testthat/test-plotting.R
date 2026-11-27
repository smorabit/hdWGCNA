is_gg <- function(x) inherits(x, 'ggplot') || inherits(x, 'patchwork')

test_that('PlotDendrogram draws without error', {
  local_null_device()
  expect_no_error(PlotDendrogram(sim_wgcna('network'), main = 'test'))
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
