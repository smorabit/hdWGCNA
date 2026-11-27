test_that('FindMajorIsoforms finds the dominant isoform of each gene', {
  so <- hdw_sim
  X <- SeuratObject::LayerData(so, layer = 'counts')
  X1 <- round(X * 0.9)
  X2 <- X - X1
  rownames(X1) <- paste0(rownames(X), '.1')
  rownames(X2) <- paste0(rownames(X), '.2')
  so[['iso']] <- SeuratObject::CreateAssay5Object(counts = rbind(X1, X2))
  SeuratObject::DefaultAssay(so) <- 'iso'
  so <- suppressWarnings(SetupForWGCNA(so, wgcna_name = 'iso', features = c(rownames(X1), rownames(X2))))

  res <- suppressMessages(FindMajorIsoforms(so, group.by = 'cell_type', replicate_col = 'Sample', assay = 'iso'))
  expect_setequal(names(res), c('A', 'B', 'C'))
  # isoform .1 carries 90% of the counts, so it is the major isoform
  expect_true(all(grepl('[.]1$', unlist(res))))

  # slot is a deprecated alias of layer
  expect_warning(
    res2 <- suppressMessages(FindMajorIsoforms(so, group.by = 'cell_type', replicate_col = 'Sample',
                                               assay = 'iso', slot = 'counts')),
    'deprecated'
  )
  expect_identical(res, res2)
})
