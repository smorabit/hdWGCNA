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

test_that('EnrichrBarPlot returns one plot per module; outdir is deprecated', {
  so <- sim_wgcna('kmes')
  so <- SetEnrichrTable(so, fake_enrichr_table(so))
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  plots <- EnrichrBarPlot(so, n_terms = 3)
  expect_named(plots, mods)
  expect_s3_class(plots[[1]], 'patchwork')
  # one panel per database, each with the top 3 terms
  expect_length(plots[[1]]$patches$plots, 1)
  expect_equal(nrow(plots[[1]][[2]]$data), 3)

  outdir <- withr::local_tempdir()
  expect_warning(EnrichrBarPlot(so, outdir = outdir, n_terms = 3), 'deprecated')
  expect_setequal(list.files(outdir), paste0(mods, '.pdf'))
})

test_that('EnrichrDotPlot returns a ggplot of the top terms', {
  so <- sim_wgcna('kmes')
  so <- SetEnrichrTable(so, fake_enrichr_table(so))
  p <- EnrichrDotPlot(so, database = 'GO_Biological_Process_2021', n_terms = 2)
  expect_s3_class(p, 'ggplot')
  expect_lte(length(unique(p$data$Term)), 2 * length(setdiff(levels(GetModules(so)$module), 'grey')))
})
