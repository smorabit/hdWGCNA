# Shared test fixtures built on the simulated dataset `hdw_sim`.
#
# The full pipeline is run once per test session and cached, so individual tests can start
# from any stage without rebuilding the network. Fixtures never write outside tempdir().

.fixture_cache <- new.env(parent = emptyenv())

# Preprocessed hdw_sim: normalized, variable features, scaled, PCA.
sim_preprocessed <- function() {
  if (is.null(.fixture_cache$pre)) {
    so <- hdw_sim
    so <- Seurat::NormalizeData(so, verbose = FALSE)
    so <- Seurat::FindVariableFeatures(so, verbose = FALSE)
    so <- Seurat::ScaleData(so, verbose = FALSE)
    so <- Seurat::RunPCA(so, npcs = 20, verbose = FALSE)
    .fixture_cache$pre <- so
  }
  .fixture_cache$pre
}

# hdWGCNA experiment 'sim' up to (and including) the requested stage, for cell type A.
# stages: 'setup', 'metacells', 'datexpr', 'powers', 'network', 'mes', 'kmes'
sim_wgcna <- function(stage = 'kmes') {
  stages <- c('setup', 'metacells', 'datexpr', 'powers', 'network', 'mes', 'kmes')
  stage <- match.arg(stage, stages)
  if (is.null(.fixture_cache$stages)) {
    .fixture_cache$stages <- list()
    # metacell sampling is random: fix the seed so the fixture does not depend on which
    # tests ran before (restore the caller's RNG state afterwards)
    withr::local_seed(12345)
    so <- sim_preprocessed()
    so <- SetupForWGCNA(so, gene_select = 'fraction', fraction = 0.05, wgcna_name = 'sim')
    .fixture_cache$stages$setup <- so
    so <- MetacellsByGroups(
      so,
      group.by = c('cell_type', 'Sample'), reduction = 'pca', k = 20, max_shared = 10,
      min_cells = 40, ident.group = 'cell_type', verbose = FALSE
    )
    so <- NormalizeMetacells(so)
    .fixture_cache$stages$metacells <- so
    so <- SetDatExpr(so, group_name = 'A', group.by = 'cell_type', assay = 'RNA', layer = 'data')
    .fixture_cache$stages$datexpr <- so
    so <- suppressWarnings(TestSoftPowers(so))
    .fixture_cache$stages$powers <- so
    invisible(utils::capture.output(
      so <- ConstructNetwork(
        so,
        tom_outdir = file.path(tempdir(), 'hdw_test_TOM'), overwrite_tom = TRUE,
        minModuleSize = 20, verbose = 0
      )
    ))
    .fixture_cache$stages$network <- so
    so <- ModuleEigengenes(so, group.by.vars = 'Sample', verbose = FALSE)
    .fixture_cache$stages$mes <- so
    so <- ModuleConnectivity(so, group.by = 'cell_type', group_name = 'A')
    .fixture_cache$stages$kmes <- so
  }
  .fixture_cache$stages[[stage]]
}

# Best-matching inferred module for each planted module of cell type A (M1, M2, M4).
# recall: share of planted genes in the matched module.
# top_precision: share of the matched module's top genes by kME (as many as were planted)
#   that are planted genes. WGCNA also attaches some uncorrelated (kME ~ 0) genes to modules,
#   so plain Jaccard varies with metacell sampling; recall and top_precision do not.
sim_module_recovery <- function(so) {
  mods <- GetModules(so)
  mods <- mods[mods$module != 'grey', ]
  truth <- hdw_sim_modules[hdw_sim_modules$module %in% c('M1', 'M2', 'M4'), ]
  do.call(rbind, lapply(split(truth$gene_name, truth$module), function(true_genes) {
    jac <- vapply(split(mods$gene_name, as.character(mods$module)), function(g) {
      length(intersect(g, true_genes)) / length(union(g, true_genes))
    }, numeric(1))
    best <- names(which.max(jac))
    cur <- mods[mods$module == best, ]
    kme_col <- paste0('kME_', best)
    top_precision <- NA_real_
    if (kme_col %in% colnames(cur)) {  # kMEs exist only after ModuleConnectivity
      top <- cur$gene_name[order(-cur[[kme_col]])][seq_along(true_genes)]
      top_precision <- mean(top %in% true_genes, na.rm = TRUE)
    }
    data.frame(
      best_match = best,
      jaccard = max(jac),
      recall = mean(true_genes %in% cur$gene_name),
      top_precision = top_precision
    )
  }))
}

# The 'kmes' stage plus a UMAP (needed by some plotting functions).
sim_with_umap <- function() {
  if (is.null(.fixture_cache$umap)) {
    withr::local_seed(12345)
    so <- sim_wgcna('kmes')
    so <- suppressWarnings(Seurat::RunUMAP(so, dims = 1:10, verbose = FALSE))
    .fixture_cache$umap <- so
  }
  .fixture_cache$umap
}

# Open a null graphics device for the rest of the calling test.
local_null_device <- function(env = parent.frame()) {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off(), envir = env)
}

# Module preservation of the cell type A network (reference) in cell type B (query).
sim_preservation <- function() {
  if (is.null(.fixture_cache$preservation)) {
    withr::local_seed(12345)
    # WGCNA::modulePreservation saves permutedStats-*.RData to the working directory
    withr::local_dir(withr::local_tempdir())
    ref <- sim_wgcna('kmes')
    query <- suppressWarnings(SetDatExpr(ref, group_name = 'B', group.by = 'cell_type', assay = 'RNA'))
    invisible(utils::capture.output(query <- suppressWarnings(
      ModulePreservation(query, seurat_ref = ref, name = 'A_in_B', n_permutations = 10, verbose = 0)
    )))
    .fixture_cache$preservation <- query
  }
  .fixture_cache$preservation
}

# A TF network on a synthetic motif matrix: three hub genes of the planted module M2 and three
# unrelated genes act as "TFs" with motifs in 60 target genes (30 M2 genes, 30 others).
sim_tf_network <- function() {
  if (is.null(.fixture_cache$tf)) {
    withr::local_seed(12345)
    so <- sim_wgcna('kmes')
    m2 <- hdw_sim_modules$gene_name[hdw_sim_modules$module == 'M2']  # sorted by loading
    noise <- setdiff(GetWGCNAGenes(so), hdw_sim_modules$gene_name)
    tfs <- c(m2[1:3], noise[1:3])
    genes <- unique(c(tfs, m2[4:33], noise[4:33]))
    motif_matrix <- matrix(
      TRUE, nrow = length(genes), ncol = length(tfs),
      dimnames = list(genes, paste0('MOTIF_', tfs))
    )
    so <- SetMotifMatrix(so, motif_matrix)
    so <- SetMotifs(so, data.frame(motif_ID = colnames(motif_matrix), motif_name = tfs, gene_name = tfs))
    model_params <- list(objective = 'reg:squarederror', max_depth = 1, eta = 0.1, nthread = 1, alpha = 0.5)
    invisible(utils::capture.output(so <- ConstructTFNetwork(so, model_params = model_params, nfold = 3)))
    .fixture_cache$tf <- list(so = so, m2_tfs = m2[1:3], noise_tfs = noise[1:3], m2 = m2)
  }
  .fixture_cache$tf
}


# Run code that calls a deprecated hdWGCNA function, muting only the deprecation warning.
suppress_deprecated <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl('is deprecated', conditionMessage(w))) invokeRestart('muffleWarning')
  })
}
