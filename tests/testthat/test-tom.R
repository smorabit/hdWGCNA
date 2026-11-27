test_that('GetTOM finds the TOM when tom_outdir is an absolute path', {
  skip_on_cran()
  so <- sim_wgcna('powers')
  tom_dir <- file.path(tempdir(), 'abs_tom_test')
  invisible(utils::capture.output(
    so <- ConstructNetwork(so, tom_outdir = tom_dir, overwrite_tom = TRUE,
                           minModuleSize = 20, verbose = 0)
  ))
  tom_file <- GetNetworkData(so)$TOMFiles
  expect_true(file.exists(tom_file))

  TOM <- GetTOM(so)
  genes <- GetModules(so)$gene_name
  expect_equal(dim(TOM), c(length(genes), length(genes)))
  expect_equal(rownames(TOM), genes)
})

test_that('GetTOM repairs TOM paths stored by older versions with an absolute tom_outdir', {
  so <- sim_wgcna('network')
  tom_file <- GetNetworkData(so)$TOMFiles

  # older versions stored paste0(getwd(), '/', <absolute path>)
  net <- GetNetworkData(so)
  net$TOMFiles <- paste0('/some/old/workdir/', tom_file)
  so <- SetNetworkData(so, net)

  expect_equal(dim(GetTOM(so))[1], nrow(GetModules(so)))
})
