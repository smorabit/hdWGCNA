test_that('ModuleTraitCorrelation works with one and with several traits', {
  so <- sim_wgcna('mes')
  so1 <- ModuleTraitCorrelation(so, traits = 'nCount_RNA', group.by = 'cell_type')
  mt1 <- GetModuleTraitCorrelation(so1)
  expect_true(all(c('cor', 'pval', 'fdr') %in% names(mt1)))
  expect_equal(rownames(mt1$cor$all_cells), 'nCount_RNA')

  so2 <- ModuleTraitCorrelation(so, traits = c('nCount_RNA', 'nFeature_RNA'), group.by = 'cell_type')
  mt2 <- GetModuleTraitCorrelation(so2)
  expect_equal(rownames(mt2$cor$all_cells), c('nCount_RNA', 'nFeature_RNA'))
  expect_true(all(c('A', 'B', 'C') %in% names(mt2$cor)))
  expect_s3_class(PlotModuleTraitCorrelation(so2), 'ggplot')
})
