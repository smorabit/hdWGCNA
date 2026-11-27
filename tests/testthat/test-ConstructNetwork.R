# check that the network information gets added to the seurat object
test_that('Network info gets put into the Seurat object', {
  data(test_seurat)

  test_seurat <- ConstructNetwork(
    test_seurat,
    soft_power = 5,
    verbose = FALSE
  )


  # # remove this in a sec:
  # test_seurat <- ModuleEigengenes(test_seurat)
  # test_seurat <- ModuleConnectivity(test_seurat)
  # head(GetModules(test_seurat))

  net <- GetNetworkData(test_seurat)

  # remove the TOM directory
  unlink('TOM', recursive = TRUE)
  # file.remove('ConsensusTOM-block.1.rda')

  expect_equal(is.null(net), FALSE)
})

# check that the TOM gets made
test_that('The TOM gets made', {
  data(test_seurat)

  wgcna_name <- 'test'

  test_seurat <- ConstructNetwork(
    test_seurat,
    soft_power = 5,
    verbose = FALSE
  )

  check <- file.exists(paste0('TOM/', wgcna_name, '_TOM.rda'))

  # remove the TOM directory
  unlink('TOM', recursive = TRUE)
  # file.remove('ConsensusTOM-block.1.rda')

  expect_equal(check, TRUE)
})

# check that SetTOM / GetTOM round-trip works with a manually supplied matrix
test_that('SetTOM stores TOM in Seurat object and GetTOM retrieves it', {
  data(test_seurat)

  # build a small named dummy matrix
  genes <- GetWGCNAGenes(test_seurat)
  n <- length(genes)
  dummy_tom <- matrix(runif(n * n), nrow = n, ncol = n)
  rownames(dummy_tom) <- genes
  colnames(dummy_tom) <- genes

  test_seurat <- SetTOM(test_seurat, dummy_tom)
  retrieved <- GetTOM(test_seurat)

  expect_equal(retrieved, dummy_tom)
})

# check that GetTOM returns the in-object TOM without needing a file on disk
test_that('GetTOM returns in-object TOM without a file on disk', {
  data(test_seurat)

  genes <- GetWGCNAGenes(test_seurat)
  n <- length(genes)
  dummy_tom <- matrix(runif(n * n), nrow = n, ncol = n)
  rownames(dummy_tom) <- genes
  colnames(dummy_tom) <- genes

  test_seurat <- SetTOM(test_seurat, dummy_tom)

  # GetTOM should return from the object without touching any file
  expect_equal(is.null(GetTOM(test_seurat)), FALSE)
})

# check that store_tom_in_seurat=TRUE stores a correctly-shaped TOM in the object
test_that('ConstructNetwork with store_tom_in_seurat=TRUE stores TOM in Seurat object', {
  data(test_seurat)

  test_seurat <- ConstructNetwork(
    test_seurat,
    soft_power = 5,
    store_tom_in_seurat = TRUE,
    verbose = FALSE
  )

  tom <- GetTOM(test_seurat)
  modules <- GetModules(test_seurat)

  # remove the TOM directory
  unlink('TOM', recursive = TRUE)

  # TOM should be a square matrix with genes matching the modules table
  expect_equal(is.matrix(tom), TRUE)
  expect_equal(nrow(tom), nrow(modules))
  expect_equal(ncol(tom), nrow(modules))
  expect_equal(rownames(tom), modules$gene_name)
  expect_equal(colnames(tom), modules$gene_name)
})

# check that store_tom_in_seurat=TRUE with saveConsensusTOMs=FALSE throws an error
test_that('store_tom_in_seurat=TRUE with saveConsensusTOMs=FALSE throws an error', {
  data(test_seurat)

  expect_error(
    ConstructNetwork(
      test_seurat,
      soft_power = 5,
      store_tom_in_seurat = TRUE,
      saveConsensusTOMs = FALSE,
      verbose = FALSE
    )
  )
})

# regression test for #397 / #341: ConstructNetwork took the modules of the active experiment
test_that('ConstructNetwork uses wgcna_name, not the active experiment', {
  so <- sim_wgcna('powers')
  # make a second experiment the active one
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

# regression test for #182: WGCNA temp files are kept in tom_outdir, and results don't change
test_that('ConstructNetwork leaves no files in the working directory and keeps results', {
  so <- sim_wgcna('powers')
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  tom_dir <- file.path(tempdir(), 'tom_182')
  invisible(utils::capture.output(
    so2 <- ConstructNetwork(so, tom_outdir = tom_dir, overwrite_tom = TRUE,
                            minModuleSize = 20, verbose = 0)
  ))
  expect_length(list.files(wd, all.files = TRUE, no.. = TRUE, recursive = TRUE), 0)
  expect_true(file.exists(file.path(tom_dir, 'sim_TOM.rda')))

  ref <- GetModules(sim_wgcna('network'))
  new <- GetModules(so2)
  expect_equal(new$gene_name, ref$gene_name)
  expect_equal(as.character(new$module), as.character(ref$module))
})

# regression test for #76: numericLabels = TRUE gave numeric module colours
test_that('ConstructNetwork converts numeric module labels to colours', {
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
