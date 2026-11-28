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

test_that('ME loadings are numeric and do not accumulate across runs', {
  skip_on_cran()
  so <- sim_wgcna('mes')
  loadings <- GetMELoadings(so, harmonized = FALSE)
  expect_true(is.numeric(loadings))
  expect_true(all(names(loadings) %in% rownames(so)))
  expect_true(is.numeric(GetMELoadings(so, harmonized = TRUE)))

  so <- ModuleEigengenes(so, group.by.vars = 'Sample', verbose = FALSE)
  expect_length(GetMELoadings(so, harmonized = FALSE), length(loadings))
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

test_that('ModuleExprScore and AvgModuleExpr score every module', {
  so <- sim_wgcna('kmes')
  mods <- setdiff(levels(GetModules(so)$module), 'grey')

  # ctrl = 20: AddModuleScore's default of 100 control genes per bin is too many for 1,000 genes
  so <- ModuleExprScore(so, n_genes = 25, method = 'Seurat', ctrl = 20)
  scores <- GetModuleScores(so)
  expect_equal(dim(scores), c(ncol(so), length(mods)))

  if (requireNamespace('UCell', quietly = TRUE)) {
    so <- ModuleExprScore(so, n_genes = 25, method = 'UCell')
    expect_equal(dim(GetModuleScores(so)), c(ncol(so), length(mods)))
  }

  so <- AvgModuleExpr(so, n_genes = 25)
  expect_equal(ncol(GetAvgModuleExpr(so)), length(mods))
})
