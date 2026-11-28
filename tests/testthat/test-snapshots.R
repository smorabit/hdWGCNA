# Visual regression tests with vdiffr. The plots are built from computed data
# (PCA, harmony, network layouts) whose last digits can differ between operating
# systems and BLAS libraries, so these run locally only (skipped on CI and CRAN).
# After an intended visual change, review and accept with testthat::snapshot_review().
expect_plot_snapshot <- function(title, p) {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed('vdiffr')
  withr::local_seed(42)
  vdiffr::expect_doppelganger(title, p)
}

test_that('dendrogram snapshot', {
  expect_plot_snapshot('PlotDendrogram', PlotDendrogram(sim_wgcna('network'), main = 'Dendrogram'))
})

test_that('radar snapshot', {
  expect_plot_snapshot('ModuleRadarPlot', ModuleRadarPlot(sim_wgcna('mes'), group.by = 'Sample', ncol = 2))
})

test_that('hub gene plot snapshots', {
  so <- sim_wgcna('kmes')
  expect_plot_snapshot('PlotKMEs', PlotKMEs(so, n_hubs = 5))
  expect_plot_snapshot('ModuleNetworkPlot', ModuleNetworkPlot(so, n_inner = 5, n_outer = 8)[[1]])
  withr::local_seed(1)
  expect_plot_snapshot('HubGeneNetworkPlot', HubGeneNetworkPlot(so, n_hubs = 4, n_other = 4))
})

test_that('module correlation network snapshot', {
  withr::local_seed(1)
  expect_plot_snapshot('ModuleCorrNetwork', ModuleCorrNetwork(sim_with_umap(), cor_cutoff = 0, label_vertices = TRUE))
})

test_that('soft power snapshot', {
  plots <- PlotSoftPowers(sim_wgcna('powers'))
  expect_plot_snapshot('PlotSoftPowers', patchwork::wrap_plots(plots, ncol = 2))
})

test_that('Enrichr bar plot snapshot', {
  so <- sim_wgcna('kmes')
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  enrich_df <- data.frame(
    Term = paste('GO term', 1:6),
    Overlap = '5/100',
    P.value = 0.001,
    Adjusted.P.value = 0.01,
    Combined.Score = c(80, 60, 45, 30, 20, 10),
    Genes = 'Gene0001;Gene0002',
    db = 'GO_Biological_Process_2021',
    module = mods[1]
  )
  so <- SetEnrichrTable(so, enrich_df)
  expect_plot_snapshot('EnrichrBarPlot', EnrichrBarPlot(so, n_terms = 5)[[1]])
})
