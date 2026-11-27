test_that('AggregatePseudobulk sums counts per replicate and group', {
  so <- hdw_sim
  X <- SeuratObject::LayerData(so, layer = 'counts')
  se <- AggregatePseudobulk(X, so@meta.data, replicate_col = 'Sample', group_col = 'cell_type')

  expect_s4_class(se, 'SummarizedExperiment')
  expect_equal(ncol(se), 24) # 8 samples x 3 cell types
  cells <- colnames(so)[so$Sample == 'S01' & so$cell_type == 'A']
  pb <- SummarizedExperiment::assay(se)
  expect_equal(
    unname(pb[, 'S01.A']),
    unname(Matrix::rowSums(X[rownames(pb), cells]))
  )
  expect_equal(SummarizedExperiment::colData(se)['S01.A', 'nCells'], length(cells))
})

test_that('ConstructPseudobulk is deprecated but still returns samples x genes', {
  so <- SetupForWGCNA(hdw_sim, gene_select = 'fraction', fraction = 0.05, wgcna_name = 'pb')
  expect_warning(
    expect_warning(
      pb <- ConstructPseudobulk(so, group.by = 'cell_type', replicate_col = 'Sample', min_reps = 2),
      'deprecated'
    ),
    'recommend at least 20 replicates'
  )
  expect_equal(nrow(pb), 24)
  expect_true(all(grepl('^[ABC]:S0[1-8]$', rownames(pb))))
  expect_true(all(colnames(pb) %in% GetWGCNAGenes(so, 'pb')))
})

test_that('pseudobulk_by_group drops only replicates with too few cells', {
  so <- hdw_sim
  drop <- colnames(so)[so$cell_type == 'C' & so$Sample == 'S03'][-(1:2)]
  so <- subset(so, cells = setdiff(colnames(so), drop))
  X <- SeuratObject::LayerData(so, layer = 'counts')

  pb <- pseudobulk_by_group(X, so@meta.data, replicate_col = 'Sample', group_col = 'cell_type')
  expect_equal(sort(names(pb)), c('A', 'B', 'C'))
  expect_equal(ncol(pb$C), 7)
  expect_false('S03' %in% colnames(pb$C))
})
