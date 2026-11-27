test_that('SetupForWGCNA keeps a user-supplied gene list', {
  genes <- rownames(hdw_sim)[1:200]
  so <- suppressWarnings(SetupForWGCNA(hdw_sim, wgcna_name = 'test', features = genes))
  expect_equal(GetWGCNAGenes(so, 'test'), genes)
  expect_equal(so@misc$active_wgcna, 'test')
})

test_that('SetupForWGCNA collapses duplicated genes', {
  genes <- rownames(hdw_sim)[1:200]
  so <- suppressWarnings(SetupForWGCNA(hdw_sim, wgcna_name = 'test', features = c(genes, genes)))
  expect_length(GetWGCNAGenes(so, 'test'), length(genes))
})

test_that('SetupForWGCNA errors for genes that are not in the object', {
  expect_error(SetupForWGCNA(hdw_sim, wgcna_name = 'test', features = paste0('Fake', 1:100)))
})

test_that('SetupForWGCNA selects genes by expression fraction or variable features', {
  so <- sim_preprocessed()
  so <- SetupForWGCNA(so, wgcna_name = 'var', gene_select = 'variable')
  expect_true(all(GetWGCNAGenes(so, 'var') %in% Seurat::VariableFeatures(so)))

  so <- SetupForWGCNA(so, wgcna_name = 'frac', gene_select = 'fraction', fraction = 0.2)
  genes <- GetWGCNAGenes(so, 'frac')
  counts <- SeuratObject::LayerData(so, layer = 'counts')
  expect_true(all(Matrix::rowMeans(counts[genes, ] > 0) >= 0.2))
  # both experiments are kept
  expect_true(all(c('var', 'frac') %in% names(so@misc)))
})
