# Opt-in tests: these hit the Enrichr web service or need large Bioconductor
# annotation packages. See helper-online.R for how to run them.

test_that('RunEnrichr finds the expected terms through the live Enrichr API', {
  skip_if_no_online_tests()
  skip_if_not_installed('enrichR')

  so <- sim_wgcna_real_genes()
  mods <- setdiff(levels(GetModules(so)$module), 'grey')[1:2]
  dbs <- 'GO_Biological_Process_2021'
  so <- suppressMessages(RunEnrichr(so, dbs = dbs, max_genes = 30, wait_time = 2))

  enrich_df <- GetEnrichrTable(so)
  expect_s3_class(enrich_df, 'data.frame')
  expect_true(all(c('Term', 'Adjusted.P.value', 'Genes', 'db', 'module') %in% colnames(enrich_df)))
  expect_setequal(unique(enrich_df$db), dbs)

  top_term <- function(mod) {
    cur <- enrich_df[enrich_df$module == mod, ]
    cur$Term[which.min(cur$Adjusted.P.value)]
  }
  expect_match(top_term(mods[1]), 'ribosom|translation|peptide', ignore.case = TRUE)
  cell_cycle_terms <- 'cell cycle|mito|chromosome|DNA replication|microtubule'
  expect_match(top_term(mods[2]), cell_cycle_terms, ignore.case = TRUE)

  # the plotting functions accept the real table
  expect_s3_class(EnrichrDotPlot(so, database = dbs, n_terms = 3), 'ggplot')
})

test_that('RunEnrichrRegulons runs positive and negative targets through Enrichr', {
  skip_if_no_online_tests()
  skip_if_not_installed('enrichR')

  so <- AssignTFRegulons(sim_tf_network()$so, strategy = 'A', reg_thresh = 0.01, n_tfs = 10)
  tfs <- unique(GetTFRegulons(so)$tf)[1:2]
  so <- SetTFRegulons(so, subset(GetTFRegulons(so), tf %in% tfs))

  # give every TF the same real target genes: ribosomal genes positively
  # correlated, cell cycle genes negatively correlated
  local_mocked_bindings(GetTFTargetGenes = function(...) {
    data.frame(
      gene = c(online_gene_sets$ribosome, online_gene_sets$cell_cycle),
      Cor = rep(c(0.5, -0.5), each = 30)
    )
  })
  so <- suppressMessages(RunEnrichrRegulons(
    so,
    dbs = 'GO_Biological_Process_2021',
    wait_time = 2
  ))

  enrich_df <- GetEnrichrRegulonTable(so)
  expect_setequal(unique(enrich_df$tf), tfs)
  expect_setequal(unique(enrich_df$target_type), c('positive', 'negative'))
  pos <- enrich_df[enrich_df$target_type == 'positive', ]
  top_pos <- pos$Term[which.min(pos$Adjusted.P.value)]
  expect_match(top_pos, 'ribosom|translation|peptide', ignore.case = TRUE)
})

test_that('MotifScan finds JASPAR motif matches in real promoters', {
  skip_if_no_online_tests()
  pkgs <- c(
    'motifmatchr', 'TFBSTools', 'JASPAR2020',
    'EnsDb.Hsapiens.v86', 'BSgenome.Hsapiens.UCSC.hg38'
  )
  for (pkg in pkgs) {
    skip_if_not_installed(pkg)
  }

  # a small Seurat object with real gene symbols
  tf_genes <- c('SPI1', 'CEBPA', 'GATA1', 'TAL1')
  genes <- unique(c(online_gene_sets$ribosome, online_gene_sets$cell_cycle, tf_genes))
  counts <- withr::with_seed(1, matrix(
    stats::rpois(length(genes) * 50, 5),
    nrow = length(genes),
    dimnames = list(genes, paste0('cell', 1:50))
  ))
  so <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  so <- SetupForWGCNA(so, wgcna_name = 'motifs', features = genes)

  pfm <- TFBSTools::getMatrixSet(
    JASPAR2020::JASPAR2020,
    opts = list(
      collection = 'CORE',
      tax_group = 'vertebrates',
      name = c('SPI1', 'CEBPA', 'GATA1', 'TAL1::TCF3'),
      all_versions = FALSE
    )
  )
  so <- suppressMessages(MotifScan(
    so,
    species_genome = 'hg38',
    pfm = pfm,
    EnsDb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  ))

  motif_mat <- GetMotifMatrix(so)
  expect_true(all(rownames(motif_mat) %in% genes))
  expect_setequal(colnames(motif_mat), names(pfm))
  expect_gt(sum(motif_mat), 0)

  motif_df <- GetMotifs(so)
  expect_true(all(c('motif_ID', 'motif_name', 'n_targets', 'gene_name') %in% colnames(motif_df)))
  # the TAL1::TCF3 dimer motif is split into its two genes, only TAL1 is in the object
  expect_true('TAL1' %in% motif_df$gene_name)
  expect_gte(length(GetMotifTargets(so)), length(pfm))
})
