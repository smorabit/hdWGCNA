test_that('BinPseudotime and PlotModuleTrajectory', {
  so <- sim_wgcna('mes')
  withr::local_seed(1)
  so$pseudotime <- stats::runif(ncol(so))
  so <- BinPseudotime(so, n_bins = 10)
  bin_cols <- grep('bin', colnames(so@meta.data), value = TRUE)
  expect_gt(length(bin_cols), 0)
  expect_equal(length(unique(so@meta.data[[bin_cols[1]]])), 10)

  p <- PlotModuleTrajectory(so, n_bins = 10)
  expect_true(inherits(p, 'ggplot') || inherits(p, 'patchwork'))
})
