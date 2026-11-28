test_that('MetacellsByGroups returns a metacell Seurat object with the expected contents', {
  so <- sim_wgcna('metacells')
  m_obj <- GetMetacellObject(so)

  expect_s4_class(m_obj, 'Seurat')
  expect_equal(rownames(m_obj), rownames(so))
  expect_true(all(c('cell_type', 'Sample') %in% colnames(m_obj@meta.data)))
  expect_true(all(as.character(Idents(m_obj)) %in% c('A', 'B', 'C')))
  # one set of metacells per cell type x sample
  expect_equal(nrow(unique(m_obj@meta.data[, c('cell_type', 'Sample')])), 24)
})

test_that('NormalizeMetacells log-normalizes the metacell counts', {
  m_obj <- GetMetacellObject(sim_wgcna('metacells'))
  expr <- SeuratObject::LayerData(m_obj, layer = 'data')
  counts <- SeuratObject::LayerData(m_obj, layer = 'counts')
  expect_false(isTRUE(all.equal(ceiling(expr), expr)))
  expect_equal(dim(expr), dim(counts))
})

# regression test for #416 / #290: a single failed group used to crash MetacellsByGroups
test_that('MetacellsByGroups removes a single group that is too small for metacells', {
  skip_on_cran()
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
  skip_on_cran()
  withr::local_options(Seurat.object.assay.version = 'v3')
  so <- SetupForWGCNA(sim_preprocessed(), gene_select = 'fraction', fraction = 0.05, wgcna_name = 'v3')
  so <- MetacellsByGroups(
    so, group.by = c('cell_type', 'Sample'), reduction = 'pca', k = 20, max_shared = 10,
    min_cells = 40, ident.group = 'cell_type', verbose = FALSE
  )
  expect_s4_class(GetMetacellObject(so)[['RNA']], 'Assay')
})

# regression test for #217: subsetting graphs without dimnames failed in MetacellsByGroups
test_that('MetacellsByGroups works with graphs that have no dimnames', {
  skip_on_cran()
  so <- Seurat::FindNeighbors(sim_preprocessed(), dims = 1:10, verbose = FALSE)
  g <- so@graphs$RNA_snn
  g@Dimnames <- list(NULL, NULL)
  so@graphs$RNA_snn <- g
  so <- SetupForWGCNA(so, gene_select = 'fraction', fraction = 0.05, wgcna_name = 'g')
  so <- MetacellsByGroups(
    so, group.by = c('cell_type', 'Sample'), reduction = 'pca', k = 20, max_shared = 10,
    min_cells = 40, ident.group = 'cell_type', verbose = FALSE
  )
  expect_gt(ncol(GetMetacellObject(so)), 0)
  expect_true('RNA_snn' %in% names(so@graphs))  # the returned object keeps its graphs
})

test_that('MetacellsByGroups requires group.by and reduction', {
  so <- sim_wgcna('setup')
  expect_error(MetacellsByGroups(so, reduction = 'pca'), 'group.by is required')
  expect_error(MetacellsByGroups(so, group.by = 'cell_type'), 'reduction is required.*pca')
  expect_error(
    MetacellsByGroups(so, group.by = c('cell_type', 'patient'), reduction = 'pca'),
    'not found in the metadata: patient'
  )
})

test_that('metacells sum the counts of their cells by default', {
  so <- sim_wgcna('metacells')
  m_obj <- GetMetacellObject(so)
  m_counts <- SeuratObject::LayerData(m_obj, layer = 'counts')
  counts <- SeuratObject::LayerData(so, layer = 'counts')
  for (i in c(1, ncol(m_obj))) {
    cells <- strsplit(m_obj$cells_merged[i], ',')[[1]]
    expect_length(cells, 20)
    expect_equal(unname(m_counts[, i]), unname(Matrix::rowSums(counts[, cells])))
  }
})

test_that('max_shared defaults to k / 5 and ident.group to the first group.by column', {
  expect_equal(default_max_shared(50), 10)
  expect_equal(default_max_shared(25), 5)

  skip_on_cran()
  so <- sim_wgcna('setup')
  so <- so[, so$cell_type == 'A' & so$Sample %in% c('S01', 'S02')]
  expect_message(
    so <- MetacellsByGroups(so, group.by = c('cell_type', 'Sample'), reduction = 'pca', k = 25),
    'Constructed [0-9]+ metacells \\(k = 25\\) in 2 group'
  )
  stats <- GetWGCNAParams(so)$metacell_stats
  expect_true(all(stats$max_shared <= 5))
  expect_true(all(stats$n <= 100))
  expect_true(all(as.character(Idents(GetMetacellObject(so))) == 'A'))
})

test_that('metacell_summary names groups that produced few metacells', {
  stats <- data.frame(name = c('A#S01', 'B#S01', 'C#S01'), n = c(100, 4, 60))
  msg <- metacell_summary(stats, k = 50)
  expect_match(msg, 'Constructed 164 metacells \\(k = 50\\) in 3 group')
  expect_match(msg, '1 group\\(s\\) produced fewer than 10 metacells.*B / S01')
  expect_no_match(metacell_summary(stats[-2, ], k = 50), 'fewer than')
})
