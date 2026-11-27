test_that('SetDatExpr stores a metacells x genes data.frame for the selected group', {
  so <- sim_wgcna('datexpr')
  datExpr <- GetDatExpr(so)
  m_obj <- GetMetacellObject(so)

  expect_s3_class(datExpr, 'data.frame')
  expect_true(all(colnames(datExpr) %in% GetWGCNAGenes(so)))
  expect_true(all(rownames(datExpr) %in% colnames(m_obj)))
  expect_true(all(m_obj$cell_type[match(rownames(datExpr), colnames(m_obj))] == 'A'))
})

test_that('SetDatExpr can return the matrix instead of the object', {
  so <- sim_wgcna('metacells')
  datExpr <- SetDatExpr(so, group_name = 'B', group.by = 'cell_type', return_seurat = FALSE)
  expect_s3_class(datExpr, 'data.frame')
  m_obj <- GetMetacellObject(so)
  expect_true(all(m_obj$cell_type[match(rownames(datExpr), colnames(m_obj))] == 'B'))
})

test_that('SetDatExpr selects several groups and can use single cells', {
  so <- sim_wgcna('metacells')
  so <- SetDatExpr(so, group_name = c('A', 'B'), group.by = 'cell_type')
  m_obj <- GetMetacellObject(so)
  groups <- m_obj$cell_type[match(rownames(GetDatExpr(so)), colnames(m_obj))]
  expect_setequal(unique(as.character(groups)), c('A', 'B'))

  so <- SetDatExpr(so, group_name = 'C', group.by = 'cell_type', use_metacells = FALSE)
  datExpr <- GetDatExpr(so)
  expect_true(all(rownames(datExpr) %in% colnames(so)))
  expect_equal(nrow(datExpr), sum(so$cell_type == 'C'))
})
