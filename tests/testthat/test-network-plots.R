test_that('edge_colors uses the module colour for within-module edges and grey otherwise', {
  modules <- data.frame(gene_name = c('a', 'b', 'c'), color = c('red', 'red', 'blue'))
  edge_df <- data.frame(Var1 = c('a', 'a', 'c'), Var2 = c('b', 'c', 'c'))
  expect_equal(edge_colors(edge_df, modules), c('red', 'grey90', 'blue'))
})

# regression test for #402: edge colours used future_sapply(), exporting the whole TOM
test_that('HubGeneNetworkPlot returns a ggraph plot without exporting large globals', {
  so <- sim_wgcna('kmes')
  withr::local_options(future.globals.maxSize = 1)  # any future export would now fail
  withr::local_seed(1)
  p <- HubGeneNetworkPlot(so, n_hubs = 3, n_other = 5)
  expect_s3_class(p, 'ggraph')
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  nodes <- p$data
  expect_equal(sum(nodes$hub), 3 * length(mods))
  expect_setequal(nodes$label[nodes$hub], as.character(GetHubGenes(so, n_hubs = 3)$gene_name))

  g <- HubGeneNetworkPlot(so, n_hubs = 3, n_other = 5, edge_prop = 1, return_graph = TRUE)
  expect_s3_class(g, 'igraph')
  expect_false(any(igraph::which_loop(g)))
  # with all edges kept, the cutoff leaves every gene at least one edge
  expect_true(all(igraph::degree(g) > 0))
  local_null_device()
  expect_no_error(print(p))
})

test_that('ModuleUMAPPlot labels hub genes and requested genes', {
  so <- sim_wgcna('kmes')
  so <- suppressMessages(suppressWarnings(
    RunModuleUMAP(so, n_hubs = 5, n_neighbors = 15, min_dist = 0.1)
  ))
  umap_df <- GetModuleUMAP(so)
  extra <- umap_df$gene[umap_df$hub == 'other'][1]
  p <- ModuleUMAPPlot(so, label_hubs = 2, label_genes = extra)
  expect_s3_class(p, 'ggraph')
  labels <- p$data$label[p$data$label != '']
  expect_true(extra %in% labels)
  expect_true(all(as.character(GetHubGenes(so, n_hubs = 2)$gene_name) %in% labels))
  expect_equal(p$data$x[match(extra, p$data$name)], umap_df$UMAP1[match(extra, umap_df$gene)])

  expect_error(ModuleUMAPPlot(so, label_genes = 'not_a_gene'), 'not_a_gene')
  expect_s3_class(ModuleUMAPPlot(so, return_graph = TRUE), 'igraph')
  local_null_device()
  expect_no_error(print(p))
})

test_that('ModuleNetworkPlot returns one plot per module and outdir is deprecated', {
  so <- sim_wgcna('kmes')
  mods <- setdiff(levels(GetModules(so)$module), 'grey')
  plots <- ModuleNetworkPlot(so, n_inner = 5, n_outer = 5, n_conns = 20)
  expect_named(plots, mods)
  expect_s3_class(plots[[1]], 'ggraph')
  # n_conns caps the number of edges; the hub genes are the nodes
  expect_lte(nrow(ggraph::get_edges()(plots[[1]]$data)), 20)
  expect_setequal(plots[[1]]$data$name, as.character(GetHubGenes(so, n_hubs = 10, mods = mods[1])$gene_name))

  outdir <- withr::local_tempdir()
  expect_warning(
    ModuleNetworkPlot(so, outdir = outdir, n_inner = 5, n_outer = 5, n_conns = 20),
    'deprecated'
  )
  expect_setequal(list.files(outdir), paste0(mods, '.pdf'))
})

test_that('ModuleCorrNetwork returns a ggraph plot of module correlations', {
  so <- sim_with_umap()
  p <- ModuleCorrNetwork(so, cor_cutoff = 0, label_vertices = TRUE)
  expect_s3_class(p, 'ggraph')
  mods <- setdiff(colnames(GetMEs(so)), 'grey')
  expect_setequal(p$data$name, mods)
  # every pair of modules is connected when cor_cutoff = 0
  n <- length(mods)
  expect_equal(nrow(ggraph::get_edges()(p$data)), n * (n - 1) / 2)
  expect_s3_class(ModuleCorrNetwork(so, reduction = NULL), 'ggraph')
  local_null_device()
  expect_no_error(print(p))
})

test_that('SavePlotList writes one file per plot or a multi-page PDF', {
  plots <- list(a = ggplot2::ggplot(), b = ggplot2::ggplot())
  outdir <- withr::local_tempdir()
  paths <- SavePlotList(plots, outdir = outdir, width = 3, height = 3)
  expect_setequal(basename(paths), c('a.pdf', 'b.pdf'))
  expect_true(all(file.exists(paths)))
  f <- file.path(outdir, 'all.pdf')
  SavePlotList(plots, file = f)
  expect_true(file.exists(f))
  expect_error(SavePlotList(plots), 'exactly one')
})
