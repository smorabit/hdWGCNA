test_that('CheckSeurat5 is TRUE for Seurat 5 and later', {
  expect_true(CheckSeurat5())
})

test_that('check_suggests passes for installed packages', {
  expect_true(check_suggests('stats'))
})

test_that('check_suggests explains how to install missing packages', {
  err <- tryCatch(
    check_suggests(c('notarealpackage1', 'notarealpackage2'), feature = 'X'),
    error = conditionMessage
  )
  expect_match(err, 'Packages .* are required for X')
  expect_match(err, "install.packages(c('notarealpackage1', 'notarealpackage2'))", fixed = TRUE)

  # Bioconductor packages get a BiocManager hint
  testthat::local_mocked_bindings(bioc_suggests = 'notarealpackage2')
  err <- tryCatch(check_suggests(c('notarealpackage1', 'notarealpackage2')), error = conditionMessage)
  expect_match(err, "BiocManager::install(c('notarealpackage2'))", fixed = TRUE)
})

test_that('check_reduction gives a clear error for a missing reduction', {
  so <- sim_preprocessed()
  expect_true(check_reduction(so, 'pca'))
  expect_error(check_reduction(so, 'umap'), 'Reduction umap not found.*pca')
})

test_that('CheckWGCNAName explains a missing experiment', {
  so <- sim_preprocessed()
  expect_error(CheckWGCNAName(so, NULL), 'Run SetupForWGCNA')
  expect_error(GetModules(so), 'Run SetupForWGCNA')
  expect_error(CheckWGCNAName(sim_wgcna('setup'), 'nope'), 'Invalid wgcna_name')
})
