# Planted modules: M1 is active in all cell types, M2 and M4 only in cell type A
# (see ?hdw_sim), so in cell type B M1 should be preserved and M2 / M4 should not.
test_that('ModulePreservation separates preserved and cell-type-specific modules', {
  skip_on_cran()
  so <- sim_preservation()
  z <- GetModulePreservation(so, 'A_in_B')$Z
  rec <- sim_module_recovery(sim_wgcna('kmes'))

  expect_gt(z[rec['M1', 'best_match'], 'Zsummary.pres'], 10)
  expect_lt(z[rec['M2', 'best_match'], 'Zsummary.pres'], 2)
  expect_lt(z[rec['M4', 'best_match'], 'Zsummary.pres'], 2)
})

test_that('preservation plots return ggplot objects', {
  skip_on_cran()
  so <- sim_preservation()
  p <- PlotModulePreservation(so, name = 'A_in_B', statistics = 'summary')
  expect_true(all(vapply(p, function(x) inherits(x, 'ggplot'), logical(1))))
  p2 <- PlotModulePreservationLollipop(so, name = 'A_in_B', features = 'Zsummary.pres', fdr = FALSE)
  expect_true(inherits(p2, 'ggplot') || all(vapply(p2, inherits, logical(1), what = 'ggplot')))
})
