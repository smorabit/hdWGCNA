test_that('ProjectModules projects modules from a reference without spurious warnings', {
  ref <- sim_wgcna('kmes')
  query <- sim_preprocessed()
  expect_no_warning(
    query <- suppressMessages(ProjectModules(query, seurat_ref = ref, wgcna_name = 'sim', wgcna_name_proj = 'proj')),
    message = 'will not be projected'
  )
  proj_mes <- GetMEs(query, harmonized = FALSE, wgcna_name = 'proj')
  ref_mods <- setdiff(levels(GetModules(ref)$module), 'grey')
  expect_true(all(ref_mods %in% colnames(proj_mes)))
  expect_equal(nrow(proj_mes), ncol(query))
})

test_that('ProjectModules accepts a user-defined modules table', {
  query <- sim_preprocessed()
  modules <- data.frame(
    gene_name = hdw_sim_modules$gene_name,
    module = factor(hdw_sim_modules$module),
    color = c(M1 = 'red', M2 = 'blue', M3 = 'green', M4 = 'orange', M5 = 'purple')[hdw_sim_modules$module]
  )
  query <- suppressMessages(suppressWarnings(
    ProjectModules(query, modules = modules, wgcna_name_proj = 'planted')
  ))
  mes <- GetMEs(query, harmonized = FALSE, wgcna_name = 'planted')
  expect_true(all(c('M1', 'M2', 'M3', 'M4', 'M5') %in% colnames(mes)))
})
