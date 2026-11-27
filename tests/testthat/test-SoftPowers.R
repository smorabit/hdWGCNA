test_that('TestSoftPowers stores a power table with the requested powers', {
  so <- sim_wgcna('datexpr')
  so <- suppressWarnings(TestSoftPowers(so, powers = c(4, 6, 8), corFnc = 'cor'))
  pt <- GetPowerTable(so)
  expect_equal(pt$Power, c(4, 6, 8))
  expect_true(all(c('SFT.R.sq', 'mean.k.') %in% colnames(pt)))
})

test_that('PlotSoftPowers returns one ggplot per statistic', {
  plots <- PlotSoftPowers(sim_wgcna('powers'))
  expect_type(plots, 'list')
  expect_true(all(vapply(plots, inherits, logical(1), what = 'ggplot')))
})
