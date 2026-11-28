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

test_that('pick_soft_power falls back to the best fit when none reaches 0.8', {
  pt <- data.frame(Power = 1:4, SFT.R.sq = c(0.2, 0.85, 0.9, 0.7))
  expect_equal(pick_soft_power(pt), 2)
  pt$SFT.R.sq <- c(0.2, 0.5, 0.6, 0.4)
  expect_message(expect_equal(pick_soft_power(pt), 3), 'best fit')
})
