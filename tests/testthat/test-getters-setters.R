# Round-trip tests for the getters and setters. Values are stored in an experiment that is
# NOT the active one, to check that wgcna_name is respected (see issues 397 / 341).

two_experiments <- function() {
  so <- SetupForWGCNA(hdw_sim, gene_select = 'fraction', fraction = 0.05, wgcna_name = 'first')
  so <- SetupForWGCNA(so, gene_select = 'fraction', fraction = 0.05, wgcna_name = 'second')
  so  # 'second' is active
}

test_that('SetupForWGCNA makes the new experiment active', {
  so <- two_experiments()
  expect_equal(GetActiveWGCNAName(so), 'second')
  so <- SetActiveWGCNA(so, 'first')
  expect_equal(GetActiveWGCNAName(so), 'first')
  expect_identical(GetActiveWGCNA(so), so@misc$first)
})

test_that('simple getters and setters round-trip and respect wgcna_name', {
  so <- two_experiments()
  value <- data.frame(a = 1:3, b = letters[1:3])
  pairs <- list(
    c('SetMetacellParams', 'GetMetacellParams'),
    c('SetWGCNAParams', 'GetWGCNAParams'),
    c('SetPowerTable', 'GetPowerTable'),
    c('SetNetworkData', 'GetNetworkData'),
    c('SetDegrees', 'GetDegrees'),
    c('SetEnrichrTable', 'GetEnrichrTable'),
    c('SetEnrichrRegulonTable', 'GetEnrichrRegulonTable'),
    c('SetModuleScores', 'GetModuleScores'),
    c('SetAvgModuleExpr', 'GetAvgModuleExpr'),
    c('SetMotifOverlap', 'GetMotifOverlap'),
    c('SetMotifScores', 'GetMotifScores'),
    c('SetModuleUMAP', 'GetModuleUMAP'),
    c('SetModuleTraitCorrelation', 'GetModuleTraitCorrelation'),
    c('SetTFNetwork', 'GetTFNetwork'),
    c('SetTFEval', 'GetTFEval'),
    c('SetTFRegulons', 'GetTFRegulons')
  )
  for (p in pairs) {
    setter <- get(p[1]); getter <- get(p[2])
    so2 <- setter(so, value, wgcna_name = 'first')
    expect_identical(getter(so2, wgcna_name = 'first'), value, label = p[2])
    # the active experiment is untouched
    expect_null(getter(so2, wgcna_name = 'second'), label = paste(p[2], 'on the active experiment'))
  }
})

test_that('getters and setters with extra arguments round-trip', {
  so <- two_experiments()
  MEs <- data.frame(blue = 1:3, grey = 3:1)

  so <- SetMEs(so, MEs, harmonized = FALSE, wgcna_name = 'first')
  so <- SetMEs(so, MEs * 2, harmonized = TRUE, wgcna_name = 'first')
  expect_identical(GetMEs(so, harmonized = FALSE, wgcna_name = 'first'), MEs)
  expect_identical(GetMEs(so, harmonized = TRUE, wgcna_name = 'first'), MEs * 2)

  so <- SetModulePreservation(so, list(x = 1), mod_name = 'query', wgcna_name = 'first')
  expect_identical(GetModulePreservation(so, mod_name = 'query', wgcna_name = 'first'), list(x = 1))

  so <- SetRegulonScores(so, MEs, target_type = 'positive', wgcna_name = 'first')
  expect_identical(GetRegulonScores(so, target_type = 'positive', wgcna_name = 'first'), MEs)
})

test_that('motif getters and setters round-trip', {
  so <- hdw_sim
  m <- matrix(TRUE, 2, 2)
  expect_identical(GetMotifMatrix(SetMotifMatrix(so, m)), m)
  df <- data.frame(motif_ID = 'x')
  expect_identical(GetMotifs(SetMotifs(so, df)), df)
  expect_identical(GetPFMList(SetPFMList(so, list(1))), list(1))
  expect_identical(GetMotifTargets(SetMotifTargets(so, list(a = 'g'))), list(a = 'g'))
})

test_that('SetModules, SetWGCNAGenes and GetHubGenes work on a real network', {
  so <- sim_wgcna('kmes')
  modules <- GetModules(so)
  so2 <- SetModules(so, modules[1:10, ])
  expect_equal(nrow(GetModules(so2)), 10)

  genes <- GetWGCNAGenes(so)[1:20]
  expect_equal(GetWGCNAGenes(SetWGCNAGenes(so, genes)), genes)

  hubs <- GetHubGenes(so, n_hubs = 5)
  mods <- setdiff(levels(modules$module), 'grey')
  expect_equal(nrow(hubs), 5 * length(mods))
  expect_setequal(as.character(unique(hubs$module)), mods)
})
