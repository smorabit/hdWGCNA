test_that('ModuleEigengenes computes MEs and harmonized MEs for every module', {
  so <- sim_wgcna('mes')
  mods <- levels(GetModules(so)$module)
  MEs <- GetMEs(so, harmonized = FALSE)
  hMEs <- GetMEs(so, harmonized = TRUE)

  expect_true(all(mods %in% colnames(MEs)))
  expect_equal(dim(hMEs), dim(MEs))
  expect_equal(rownames(MEs), colnames(so))
  expect_false(isTRUE(all.equal(hMEs, MEs)))
})

# regression test for #433: the ScaleData check failed when @commands had no ScaleData entry
test_that('harmonized MEs do not depend on ScaleData having been recorded', {
  skip_on_cran()
  so <- sim_wgcna('network')
  so@commands <- list()
  so <- ModuleEigengenes(so, group.by.vars = 'Sample', verbose = FALSE)
  hMEs <- GetMEs(so, harmonized = TRUE)
  expect_equal(nrow(hMEs), ncol(so))
  expect_true(all(levels(GetModules(so)$module) %in% colnames(hMEs)))
})
