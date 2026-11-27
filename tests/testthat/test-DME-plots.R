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
