test_that('edge_colors uses the module colour for within-module edges and grey otherwise', {
  modules <- data.frame(gene_name = c('a', 'b', 'c'), color = c('red', 'red', 'blue'))
  edge_df <- data.frame(Var1 = c('a', 'a', 'c'), Var2 = c('b', 'c', 'c'))
  expect_equal(edge_colors(edge_df, modules), c('red', 'grey90', 'blue'))
})

# regression test for #402: edge colours used future_sapply(), exporting the whole TOM
test_that('HubGeneNetworkPlot runs without exporting large globals', {
  so <- sim_wgcna('kmes')
  withr::local_options(future.globals.maxSize = 1)  # any future export would now fail
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  expect_no_error(suppressMessages(HubGeneNetworkPlot(so, n_hubs = 3, n_other = 5)))
})

# regression test for #386: vertex.label.cex was hard-coded to 0
test_that('ModuleUMAPPlot passes vertex.label.cex to igraph', {
  so <- sim_wgcna('kmes')
  so <- suppressMessages(suppressWarnings(
    RunModuleUMAP(so, n_hubs = 5, n_neighbors = 15, min_dist = 0.1)
  ))
  captured <- NULL
  testthat::local_mocked_bindings(
    plot = function(x, ...) captured <<- list(...),
    .package = 'base'
  )
  suppressMessages(ModuleUMAPPlot(so, vertex.label.cex = 0.8))
  expect_equal(captured$vertex.label.cex, 0.8)
})
