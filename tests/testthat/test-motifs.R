test_that('OverlapModulesMotifs tests TF target sets against modules', {
  skip_if_not_installed('GeneOverlap')
  so <- sim_wgcna('kmes')
  m2 <- hdw_sim_modules$gene_name[hdw_sim_modules$module == 'M2']
  noise <- setdiff(GetWGCNAGenes(so), hdw_sim_modules$gene_name)
  # TF_A targets the planted module M2, TF_B targets unrelated genes
  so <- SetMotifTargets(so, list(MOTIF_A = m2[1:40], MOTIF_B = noise[1:40]))
  so <- SetMotifs(so, data.frame(motif_ID = c('MOTIF_A', 'MOTIF_B'), motif_name = c('TF_A', 'TF_B'),
                                 gene_name = c('TF_A', 'TF_B')))
  so <- OverlapModulesMotifs(so)
  overlaps <- GetMotifOverlap(so)
  expect_true(all(c('module', 'tf', 'odds_ratio', 'pval', 'fdr') %in% colnames(overlaps)))

  m2_mod <- sim_module_recovery(so)['M2', 'best_match']
  a <- overlaps[overlaps$module == m2_mod & overlaps$tf == 'MOTIF_A', ]
  b <- overlaps[overlaps$module == m2_mod & overlaps$tf == 'MOTIF_B', ]
  expect_lt(a$fdr, 0.05)
  expect_gt(a$odds_ratio, b$odds_ratio)
})
