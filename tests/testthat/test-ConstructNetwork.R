test_that('ConstructNetwork stores the network and writes the TOM', {
  so <- sim_wgcna('network')
  net <- GetNetworkData(so)
  expect_false(is.null(net))
  expect_true(file.exists(net$TOMFiles))
  expect_equal(nrow(GetModules(so)), length(GetWGCNAGenes(so)))
})

test_that('SetTOM stores a TOM in the object and GetTOM returns it without a file', {
  so <- sim_wgcna('network')
  genes <- GetWGCNAGenes(so)
  n <- length(genes)
  dummy_tom <- matrix(stats::runif(n * n), nrow = n, ncol = n, dimnames = list(genes, genes))

  so <- SetTOM(so, dummy_tom)
  net <- GetNetworkData(so)
  net$TOMFiles <- file.path(tempdir(), 'does_not_exist.rda')
  so <- SetNetworkData(so, net)
  expect_equal(GetTOM(so), dummy_tom)
})

test_that('store_tom_in_seurat = TRUE requires saveConsensusTOMs = TRUE', {
  expect_error(
    ConstructNetwork(sim_wgcna('powers'), store_tom_in_seurat = TRUE, saveConsensusTOMs = FALSE),
    'requires saveConsensusTOMs'
  )
})

# also a regression test for #182: WGCNA temp files are kept in tom_outdir
test_that('ConstructNetwork can store the TOM in the object and leaves no files in the working directory', {
  skip_on_cran()
  so <- sim_wgcna('powers')
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  tom_dir <- file.path(tempdir(), 'tom_store')
  invisible(utils::capture.output(
    so2 <- ConstructNetwork(so, tom_outdir = tom_dir, overwrite_tom = TRUE, minModuleSize = 20,
                            store_tom_in_seurat = TRUE, verbose = 0)
  ))
  expect_length(list.files(wd, all.files = TRUE, no.. = TRUE, recursive = TRUE), 0)
  expect_true(file.exists(file.path(tom_dir, 'sim_TOM.rda')))

  # same modules as the cached network built with the same inputs
  ref <- GetModules(sim_wgcna('network'))
  new <- GetModules(so2)
  expect_equal(new$gene_name, ref$gene_name)
  expect_equal(as.character(new$module), as.character(ref$module))

  # TOM stored in the object: square, named by the module genes
  tom <- GetTOM(so2)
  expect_true(is.matrix(tom))
  expect_equal(dim(tom), c(nrow(new), nrow(new)))
  expect_equal(rownames(tom), new$gene_name)
})

# regression test for #397 / #341: ConstructNetwork took the modules of the active experiment
test_that('ConstructNetwork uses wgcna_name, not the active experiment', {
  skip_on_cran()
  so <- sim_wgcna('powers')
  so <- SetupForWGCNA(so, gene_select = 'fraction', fraction = 0.05, wgcna_name = 'other')
  expect_equal(so@misc$active_wgcna, 'other')

  invisible(utils::capture.output(
    so <- ConstructNetwork(
      so, wgcna_name = 'sim', tom_outdir = file.path(tempdir(), 'tom_397'),
      overwrite_tom = TRUE, minModuleSize = 20, verbose = 0
    )
  ))
  mods <- GetModules(so, wgcna_name = 'sim')
  expect_equal(nrow(mods), length(GetWGCNAGenes(so, 'sim')))
  expect_true(length(setdiff(unique(mods$module), 'grey')) >= 3)
})

# regression test for #76: numericLabels = TRUE gave numeric module colours
test_that('ConstructNetwork converts numeric module labels to colours', {
  skip_on_cran()
  so <- sim_wgcna('powers')
  expect_message(invisible(utils::capture.output(
    so <- ConstructNetwork(so, tom_outdir = file.path(tempdir(), 'tom_76'), overwrite_tom = TRUE,
                           minModuleSize = 20, verbose = 0, numericLabels = TRUE)
  )), 'numeric module labels')
  mods <- GetModules(so)
  expect_true('grey' %in% mods$color)
  expect_false(any(grepl('^[0-9]+$', mods$color)))
  expect_setequal(as.character(mods$module), as.character(GetModules(sim_wgcna('network'))$module))
})
