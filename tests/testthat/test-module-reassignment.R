test_that('ReassignModules moves selected genes to new modules', {
  so <- sim_wgcna('kmes')
  modules <- GetModules(so)
  mods <- setdiff(levels(modules$module), 'grey')
  genes <- modules$gene_name[modules$module == mods[1]][1:5]

  # moving genes that are already in a module requires ignore = TRUE
  expect_error(ReassignModules(so, features = genes, new_modules = rep(mods[2], 5)), 'ignore=TRUE')
  so2 <- ReassignModules(so, features = genes, new_modules = rep(mods[2], 5), ignore = TRUE)
  new <- GetModules(so2)
  expect_true(all(new$module[match(genes, new$gene_name)] == mods[2]))
  expect_true(all(new$color[match(genes, new$gene_name)] == new$color[new$module == mods[2]][1]))
  expect_equal(nrow(new), nrow(modules))
})

test_that('ReassignModules rejects unknown genes and modules', {
  so <- sim_wgcna('kmes')
  expect_error(ReassignModules(so, features = 'not_a_gene', new_modules = 'grey'), 'not found')
  gene <- GetModules(so)$gene_name[1]
  expect_error(ReassignModules(so, features = gene, new_modules = 'not_a_module'), 'invalid')
})

test_that('ResetModuleNames and ResetModuleColors rename modules', {
  so <- sim_wgcna('kmes')
  so <- ResetModuleNames(so, new_name = 'SIM')
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  expect_true(all(grepl('^SIM', mods)))
  expect_true(all(mods %in% colnames(GetMEs(so))))

  new_colors <- setNames(grDevices::rainbow(length(mods)), mods)
  so <- ResetModuleColors(so, new_colors)
  modules <- GetModules(so)
  expect_equal(unname(modules$color[match(mods, modules$module)]), unname(new_colors))
})
