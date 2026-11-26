# Simulation of single-cell data with planted co-expression modules for hdWGCNA
# testing and benchmarking. Not part of the package build (data-raw/ is .Rbuildignore'd).
#
# Requires: splatter (Bioconductor), Seurat, Matrix
#
# Approach
#   1. splatter::splatSimulate() generates cell types (groups) and samples (batches) with
#      realistic mean-variance structure, but no within-cell-type gene-gene covariance.
#   2. We take splatter's per-cell expected expression ("CellMeans"), and for each planted
#      module multiply the means of its genes by a shared per-cell latent factor
#      exp(loading_g * sigma * z_c), restricted to the cell types (and samples) where the
#      module is active. z_c ~ N(shift, 1), with an optional shift by condition.
#   3. Counts are re-drawn as Poisson(means).
#
# Ground truth returned alongside the Seurat object: gene -> module, per-gene loading
# (true "hubness"), and where each module is active.

#' Default planted-module design
#'
#' M1: global program, active in all cell types
#' M2: cell type A specific
#' M3: cell type B specific
#' M4: cell type A specific, up-regulated in "disease" (true differential module)
#' M5: cell type C specific, active only in "control" samples (not preserved in disease)
default_module_design <- function() {
  data.frame(
    module = paste0('M', 1:5),
    cell_types = c('A,B,C', 'A', 'B', 'A', 'C'),
    conditions = c(
      'control,disease', 'control,disease', 'control,disease',
      'control,disease', 'control'
    ),
    disease_shift = c(0, 0, 0, 1, 0), # shift of the latent factor mean in disease samples
    stringsAsFactors = FALSE
  )
}

simulate_hdwgcna_data <- function(
  n_genes = 1000,
  cells_per_sample = 500,
  n_samples = 8,
  cell_type_prob = c(A = 0.45, B = 0.35, C = 0.20),
  module_design = default_module_design(),
  genes_per_module = 50,
  sigma = 1.2,
  sample_effect = 0,
  lib_loc = 7.5,
  seed = 12345
) {
  stopifnot(n_samples %% 2 == 0)
  set.seed(seed)

  params <- splatter::newSplatParams(
    nGenes = n_genes,
    batchCells = rep(cells_per_sample, n_samples),
    # per-gene sample (batch) effects; 0 removes them. Even small random per-gene sample
    # effects create many spurious modules after metacell aggregation within samples.
    batch.facLoc = max(sample_effect, 1e-6),
    batch.facScale = max(sample_effect, 1e-6),
    batch.rmEffect = sample_effect == 0,
    group.prob = unname(cell_type_prob),
    de.prob = 0.08,
    de.facLoc = 0.3,
    lib.loc = lib_loc,
    seed = seed
  )
  sim <- splatter::splatSimulate(params, method = 'groups', verbose = FALSE)

  means <- SummarizedExperiment::assay(sim, 'CellMeans')
  cd <- as.data.frame(SummarizedExperiment::colData(sim))

  # metadata -------------------------------------------------------------------
  ct_names <- names(cell_type_prob)
  cell_type <- ct_names[as.integer(sub('Group', '', cd$Group))]
  sample_idx <- as.integer(sub('Batch', '', cd$Batch))
  sample <- sprintf('S%02d', sample_idx)
  condition <- ifelse(sample_idx <= n_samples / 2, 'control', 'disease')

  gene_names <- sprintf('Gene%04d', seq_len(n_genes))
  rownames(means) <- gene_names

  # choose module genes among reasonably expressed genes ----------------------
  base_mean <- SummarizedExperiment::rowData(sim)$BaseGeneMean
  candidates <- which(base_mean > stats::quantile(base_mean, 0.4))
  n_mod <- nrow(module_design)
  stopifnot(length(candidates) >= n_mod * genes_per_module)
  picked <- sample(candidates, n_mod * genes_per_module)

  truth <- data.frame(
    gene_name = gene_names[picked],
    module = rep(module_design$module, each = genes_per_module),
    loading = stats::runif(n_mod * genes_per_module, 0.5, 1),
    stringsAsFactors = FALSE
  )

  # plant modules ----------------------------------------------------------------
  latent <- matrix(0,
    nrow = ncol(means), ncol = n_mod,
    dimnames = list(colnames(means), module_design$module)
  )
  for (i in seq_len(n_mod)) {
    md <- module_design[i, ]
    active <- cell_type %in% strsplit(md$cell_types, ',')[[1]] &
      condition %in% strsplit(md$conditions, ',')[[1]]
    z <- stats::rnorm(sum(active), mean = md$disease_shift * (condition[active] == 'disease'))
    latent[active, i] <- z

    idx <- picked[truth$module == md$module]
    load <- truth$loading[truth$module == md$module]
    # log-normal multiplicative factor (mean-corrected so average expression is unchanged)
    f <- exp(outer(load * sigma, z) - (load * sigma)^2 / 2)
    means[idx, active] <- means[idx, active] * f
  }

  counts <- matrix(stats::rpois(length(means), means),
    nrow = nrow(means),
    dimnames = list(gene_names, paste0('cell', seq_len(ncol(means))))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  rownames(latent) <- colnames(counts)

  meta <- data.frame(
    cell_type = factor(cell_type, levels = ct_names),
    Sample = factor(sample),
    condition = factor(condition, levels = c('control', 'disease')),
    row.names = colnames(counts)
  )

  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts, meta.data = meta,
    project = 'hdWGCNA_sim'
  )

  list(
    seurat = seurat_obj,
    modules = truth,
    design = module_design,
    latent = latent
  )
}

# standard Seurat preprocessing used for the shipped object and the benchmarks
preprocess_sim <- function(seurat_obj, npcs = 20, verbose = FALSE) {
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = verbose)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = verbose)
  seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = verbose)
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = npcs, verbose = verbose)
  seurat_obj
}

# For each true module, the best-matching inferred module and its Jaccard index.
module_recovery <- function(inferred, truth, modules = unique(truth$module)) {
  inferred <- inferred[inferred$module != 'grey', ]
  do.call(rbind, lapply(modules, function(m) {
    true_genes <- truth$gene_name[truth$module == m]
    jac <- vapply(split(inferred$gene_name, as.character(inferred$module)), function(g) {
      length(intersect(g, true_genes)) / length(union(g, true_genes))
    }, numeric(1))
    if (length(jac) == 0) {
      return(data.frame(true_module = m, best_match = NA, jaccard = 0, precision = 0, recall = 0))
    }
    best <- names(which.max(jac))
    found <- inferred$gene_name[inferred$module == best]
    data.frame(
      true_module = m, best_match = best, jaccard = max(jac),
      precision = length(intersect(found, true_genes)) / length(found),
      recall = length(intersect(found, true_genes)) / length(true_genes)
    )
  }))
}
