# RunEnrichr needs the Enrichr web service, so these tests use a synthetic results table.
fake_enrichr_table <- function(so) {
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  dbs <- c('GO_Biological_Process_2021', 'GO_Molecular_Function_2021')
  grid <- expand.grid(module = mods, db = dbs, term_id = 1:6, stringsAsFactors = FALSE)
  withr::with_seed(1, data.frame(
    Term = paste('term', grid$term_id, grid$module, grid$db),
    Overlap = '5/100',
    P.value = stats::runif(nrow(grid), 0, 0.01),
    Adjusted.P.value = stats::runif(nrow(grid), 0, 0.04),
    Combined.Score = stats::runif(nrow(grid), 5, 100),
    Genes = 'Gene0001;Gene0002',
    db = grid$db,
    module = factor(grid$module, levels = levels(GetModules(so)$module))
  ))
}

test_that('RunEnrichr requires enrichR', {
  if (!requireNamespace('enrichR', quietly = TRUE)) {
    expect_error(RunEnrichr(sim_wgcna('kmes')), 'enrichR')
  } else {
    succeed()
  }
})

test_that('EnrichrBarPlot writes one PDF per module', {
  so <- sim_wgcna('kmes')
  so <- SetEnrichrTable(so, fake_enrichr_table(so))
  outdir <- withr::local_tempdir()
  suppressMessages(EnrichrBarPlot(so, outdir = outdir, n_terms = 3))
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  expect_setequal(list.files(outdir), paste0(mods, '.pdf'))
})

test_that('EnrichrDotPlot returns a ggplot of the top terms', {
  so <- sim_wgcna('kmes')
  so <- SetEnrichrTable(so, fake_enrichr_table(so))
  p <- EnrichrDotPlot(so, database = 'GO_Biological_Process_2021', n_terms = 2)
  expect_s3_class(p, 'ggplot')
  expect_lte(length(unique(p$data$Term)), 2 * length(setdiff(levels(GetModules(so)$module), 'grey')))
})
