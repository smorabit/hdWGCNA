dmes_disease_vs_control <- function(so) {
  cells_a <- so$cell_type == 'A'
  dmes <- FindDMEs(
    so,
    barcodes1 = colnames(so)[cells_a & so$condition == 'disease'],
    barcodes2 = colnames(so)[cells_a & so$condition == 'control']
  )
  dmes$group <- 'disease_vs_control'
  dmes
}

test_that('PlotDMEsVolcano returns a ggplot', {
  so <- sim_wgcna('mes')
  p <- PlotDMEsVolcano(so, dmes_disease_vs_control(so))
  expect_s3_class(p, 'ggplot')
})

test_that('PlotDMEsLollipop works without wgcna_name and pvalue', {
  so <- sim_wgcna('mes')
  dmes <- dmes_disease_vs_control(so)
  plots <- suppressMessages(PlotDMEsLollipop(so, dmes, group.by = 'group', comparison = 'disease_vs_control'))
  expect_type(plots, 'list')
  expect_s3_class(plots[[1]], 'ggplot')
  # renders (including the background stripes) without error
  expect_silent(ggplot2::ggplot_build(plots[[1]]))
})

test_that('ModuleCorrNetwork explains a missing reduction', {
  so <- sim_wgcna('mes')
  expect_error(ModuleCorrNetwork(so, reduction = 'umap'), 'Reduction umap not found')
})

# regression test for #424: latent.vars given as metadata column names
test_that('FindDMEs accepts latent.vars as metadata column names', {
  skip_if_not_installed('lmtest')  # needed by Seurat for test.use = 'LR'
  so <- sim_wgcna('mes')
  cells_a <- so$cell_type == 'A'
  b1 <- colnames(so)[cells_a & so$condition == 'disease']
  b2 <- colnames(so)[cells_a & so$condition == 'control']

  dmes <- FindDMEs(so, barcodes1 = b1, barcodes2 = b2, test.use = 'LR', latent.vars = 'nCount_RNA')
  expect_true(all(c('p_val', 'avg_log2FC', 'p_val_adj', 'module') %in% colnames(dmes)))
  expect_gt(nrow(dmes), 0)

  expect_error(
    FindDMEs(so, barcodes1 = b1, barcodes2 = b2, test.use = 'LR', latent.vars = 'not_a_column'),
    'latent.vars not found'
  )
})

test_that('FindAllDMEs accepts latent.vars as metadata column names', {
  skip_if_not_installed('lmtest')  # needed by Seurat for test.use = 'LR'
  so <- sim_wgcna('mes')
  dmes <- suppressMessages(FindAllDMEs(so, group.by = 'cell_type', test.use = 'LR', latent.vars = 'nCount_RNA'))
  expect_setequal(unique(dmes$group), c('A', 'B', 'C'))
})

test_that('PlotLollipop orders modules by fold change and uses the given column', {
  modules <- data.frame(module = c('a', 'b', 'c'), color = c('red', 'blue', 'green'))
  dmes <- data.frame(module = c('a', 'b', 'c'), lfc = c(0.5, -1, 2), p_val_adj = c(0.01, 0.2, 0.001))
  p <- PlotLollipop(modules, dmes, pvalue = 'p_val_adj', avg_log2FC = 'lfc')
  expect_equal(levels(p$data$module), c('b', 'a', 'c'))
  expect_equal(rlang::as_label(p$mapping$x), 'lfc')
})
