# End-to-end correctness on the simulated data: the standard pipeline should recover the
# modules planted in cell type A (see ?hdw_sim).

test_that('the pipeline recovers the planted modules of cell type A', {
  so <- sim_wgcna('kmes')
  rec <- sim_module_recovery(so)

  # each planted module (M1, M2, M4) maps to a distinct inferred module that contains all
  # of its genes, and whose top genes by kME are the planted genes
  expect_equal(length(unique(rec$best_match)), 3)
  expect_true(all(rec$recall >= 0.95))
  expect_true(all(rec$top_precision >= 0.9))
})

test_that('the hub genes of a recovered module are the planted high-loading genes', {
  so <- sim_wgcna('kmes')
  rec <- sim_module_recovery(so)
  m2 <- rec['M2', 'best_match']

  hubs <- GetHubGenes(so, n_hubs = 10)
  hubs <- hubs$gene_name[hubs$module == m2]
  truth <- hdw_sim_modules[hdw_sim_modules$module == 'M2', ]

  # hub genes should be planted genes with above-median loadings
  expect_true(all(hubs %in% truth$gene_name))
  expect_gt(mean(truth$loading[match(hubs, truth$gene_name)]), stats::median(truth$loading))
})

test_that('FindDMEs finds the module planted as up-regulated in disease', {
  so <- sim_wgcna('mes')
  rec <- sim_module_recovery(so)
  m4 <- rec['M4', 'best_match']

  cells_a <- so$cell_type == 'A'
  dmes <- FindDMEs(
    so,
    barcodes1 = colnames(so)[cells_a & so$condition == 'disease'],
    barcodes2 = colnames(so)[cells_a & so$condition == 'control']
  )
  dmes <- dmes[order(-dmes$avg_log2FC), ]

  expect_equal(dmes$module[1], m4)
  expect_lt(dmes$p_val_adj[1], 0.05)
})
