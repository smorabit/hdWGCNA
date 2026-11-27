test_that('ConstructTFNetwork gives planted TFs high importance for their module targets', {
  skip_on_cran()
  skip_if_not_installed('xgboost')
  tf <- sim_tf_network()
  tf_net <- GetTFNetwork(tf$so)
  expect_true(all(c('tf', 'gene', 'Gain', 'Cor') %in% colnames(tf_net)))

  m2_targets <- tf_net[tf_net$gene %in% tf$m2, ]
  gain_m2 <- mean(m2_targets$Gain[m2_targets$tf %in% tf$m2_tfs])
  gain_noise <- mean(m2_targets$Gain[m2_targets$tf %in% tf$noise_tfs])
  expect_gt(gain_m2, 5 * gain_noise)
})

test_that('regulons, regulon scores and differential regulons', {
  skip_on_cran()
  skip_if_not_installed('xgboost')
  skip_if_not_installed('UCell')
  tf <- sim_tf_network()
  so <- AssignTFRegulons(tf$so, strategy = 'A', reg_thresh = 0.01, n_tfs = 10)
  regulons <- GetTFRegulons(so)
  expect_true(all(tf$m2_tfs %in% regulons$tf))
  targets <- GetTFTargetGenes(so, selected_tfs = tf$m2_tfs[1], depth = 1)
  expect_true(all(targets$tf == tf$m2_tfs[1]))
  expect_true(all(targets$depth == 1))
  expect_true(all(targets$gene %in% rownames(so)))

  so <- RegulonScores(so, target_type = 'positive')
  so <- RegulonScores(so, target_type = 'negative')
  expect_equal(nrow(GetRegulonScores(so, target_type = 'positive')), ncol(so))

  # regression test: FindDifferentialRegulons referred to undefined group1 / group2
  cells_a <- so$cell_type == 'A'
  dregs <- FindDifferentialRegulons(
    so,
    barcodes1 = colnames(so)[cells_a & so$condition == 'disease'],
    barcodes2 = colnames(so)[cells_a & so$condition == 'control']
  )
  expect_true('tf' %in% colnames(dregs))
  expect_s3_class(PlotDifferentialRegulons(so, dregs), 'ggplot')
  expect_s3_class(RegulonBarPlot(so, selected_tf = tf$m2_tfs[1]), 'ggplot')
})

test_that('TF network and module regulatory plots build without attaching ggraph', {
  skip_on_cran()
  skip_if_not_installed('xgboost')
  skip_if_not_installed('tidyr')
  tf <- sim_tf_network()
  so <- AssignTFRegulons(tf$so, strategy = 'A', reg_thresh = 0.01, n_tfs = 10)
  local_null_device()

  # regression test: the edge colour guide was given by name, which needs ggraph attached
  expect_no_error(print(TFNetworkPlot(so, selected_tfs = tf$m2_tfs[1])))

  mrn <- ModuleRegulatoryNetwork(so)
  expect_s3_class(mrn, 'data.frame')
  expect_no_error(print(ModuleRegulatoryHeatmap(so)))
  expect_no_error(print(ModuleRegulatoryNetworkPlot(so, layout = 'fr')))
})
