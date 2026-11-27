# test that the output format is Seurat
test_that('Metacell output is Seurat format', {
  data(test_seurat)

  test_seurat <- SetupForWGCNA(
    test_seurat,
    wgcna_name = 'test',
    features = rownames(test_seurat)
  )

  test_seurat <- MetacellsByGroups(
    test_seurat,
    group.by = c('cell_type'),
    k = 5,
    max_shared = 5,
    ident.group = 'cell_type',
    target_metacells = 50,
    min_cells = 10
  )

  m_obj <- GetMetacellObject(test_seurat)

  expect_equal('Seurat' %in% class(m_obj), TRUE)
})

# test that the genes are the same in the metacell and the main seurat object
test_that('Metacell features are same as Seurat object', {
  data(test_seurat)

  test_seurat <- SetupForWGCNA(
    test_seurat,
    wgcna_name = 'test',
    features = rownames(test_seurat)
  )

  test_seurat <- MetacellsByGroups(
    test_seurat,
    group.by = c('cell_type'),
    k = 5,
    max_shared = 5,
    ident.group = 'cell_type',
    target_metacells = 50,
    min_cells = 10
  )

  m_obj <- GetMetacellObject(test_seurat)

  expect_equal(rownames(m_obj), rownames(test_seurat))
})


# test that the meta-data is getting passed
test_that('Meta-data gets passed to the metacell object', {
  data(test_seurat)

  test_seurat <- SetupForWGCNA(
    test_seurat,
    wgcna_name = 'test',
    features = rownames(test_seurat)
  )

  test_seurat <- MetacellsByGroups(
    test_seurat,
    group.by = c('cell_type', 'annotation'),
    k = 5,
    max_shared = 5,
    ident.group = 'cell_type',
    target_metacells = 50,
    min_cells = 10
  )

  m_obj <- GetMetacellObject(test_seurat)

  expect_equal(all(c('cell_type', 'annotation') %in% colnames(m_obj@meta.data)), TRUE)
})

# test that the idents are expected
test_that('Metacell object Idents are expected', {
  data(test_seurat)

  test_seurat <- SetupForWGCNA(
    test_seurat,
    wgcna_name = 'test',
    features = rownames(test_seurat)
  )

  test_seurat <- MetacellsByGroups(
    test_seurat,
    group.by = c('Sample', 'cell_type', 'annotation'),
    k = 5,
    max_shared = 5,
    ident.group = 'Sample',
    target_metacells = 10,
    min_cells = 10
  )

  m_obj <- GetMetacellObject(test_seurat)

  check <- all(as.character(Idents(m_obj)) %in% as.character(test_seurat$Sample))
  expect_equal(check, TRUE)
})

# test that normalizing the metacells works as intended
test_that('Normalizing the metacell object works', {
  data(test_seurat)

  test_seurat <- SetupForWGCNA(
    test_seurat,
    wgcna_name = 'test',
    features = rownames(test_seurat)
  )

  test_seurat <- MetacellsByGroups(
    test_seurat,
    group.by = c('cell_type'),
    k = 5,
    max_shared = 5,
    ident.group = 'cell_type',
    target_metacells = 50,
    min_cells = 10
  )

  test_seurat <- NormalizeMetacells(test_seurat, verbose = FALSE)
  m_obj <- GetMetacellObject(test_seurat)
  if (CheckSeurat5()) {
    expr <- SeuratObject::LayerData(test_seurat, layer = 'data')
  } else {
    expr <- Seurat::GetAssayData(test_seurat, slot = 'data')
  }

  check <- all.equal(ceiling(expr), expr) == TRUE
  expect_equal(check, FALSE)
})

# regression test for #416 / #290: a single failed group used to crash MetacellsByGroups
test_that('MetacellsByGroups removes a single group that is too small for metacells', {
  so <- sim_preprocessed()
  small_group <- so$cell_type == 'C' & so$Sample == 'S03'
  so <- subset(so, cells = c(colnames(so)[!small_group], colnames(so)[small_group][1:30]))
  so <- SetupForWGCNA(so, gene_select = 'fraction', fraction = 0.05, wgcna_name = 'mc')

  expect_warning(
    so <- suppressMessages(MetacellsByGroups(
      so, group.by = c('cell_type', 'Sample'), reduction = 'pca', k = 25, max_shared = 10,
      min_cells = 20, ident.group = 'cell_type', verbose = FALSE
    )),
    'failed for 1 group.*C#S03'
  )
  m_obj <- GetMetacellObject(so)
  expect_false(any(m_obj$cell_type == 'C' & m_obj$Sample == 'S03'))
  expect_gt(ncol(m_obj), 0)
})

test_that('MetacellsByGroups stops clearly when every group fails', {
  so <- sim_preprocessed()
  so <- SetupForWGCNA(so, gene_select = 'fraction', fraction = 0.05, wgcna_name = 'mc')
  expect_error(
    suppressWarnings(MetacellsByGroups(
      so, group.by = c('cell_type', 'Sample'), reduction = 'pca', k = 400, max_shared = 1,
      min_cells = 20, ident.group = 'cell_type', verbose = FALSE
    )),
    'failed for all groups'
  )
})

# regression test for #296: JoinLayers was called on v3 assays
test_that('MetacellsByGroups works when Seurat creates v3 assays', {
  withr::local_options(Seurat.object.assay.version = 'v3')
  so <- SetupForWGCNA(sim_preprocessed(), gene_select = 'fraction', fraction = 0.05, wgcna_name = 'v3')
  so <- MetacellsByGroups(
    so, group.by = c('cell_type', 'Sample'), reduction = 'pca', k = 20, max_shared = 10,
    min_cells = 40, ident.group = 'cell_type', verbose = FALSE
  )
  expect_s4_class(GetMetacellObject(so)[['RNA']], 'Assay')
})
