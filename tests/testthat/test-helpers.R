test_that('CheckSeurat5 is TRUE for Seurat 5 and later', {
  expect_true(CheckSeurat5())
})

test_that('check_suggests passes for installed packages', {
  expect_true(check_suggests('stats'))
})

test_that('check_suggests explains how to install missing packages', {
  err <- tryCatch(
    check_suggests(c('notarealpackage1', 'GenomicRanges_not_real', 'ggradar'), feature = 'X'),
    error = conditionMessage
  )
  expect_match(err, 'required for X')
  expect_match(err, "install.packages(c('notarealpackage1'", fixed = TRUE)
  if (!requireNamespace('ggradar', quietly = TRUE)) {
    expect_match(err, "remotes::install_github('ricardo-bion/ggradar')", fixed = TRUE)
  }
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
