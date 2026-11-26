# Generates the simulated example dataset shipped with hdWGCNA:
#   data/hdw_sim.rda          Seurat object (raw counts + metadata)
#   data/hdw_sim_modules.rda  ground-truth planted modules
#
# Run from the package root:  pixi run Rscript data-raw/hdw_sim.R
#
# Default parameters were chosen (2026-11-26) so that the standard pipeline
# (metacells by cell_type x Sample, k = 20; signed network) recovers the planted modules
# of cell type A with recall 1.0 and Jaccard ~0.9, while the data stays ~1.2 MB.

source("data-raw/sim_functions.R")

sim <- simulate_hdwgcna_data(
  n_genes = 1000,
  cells_per_sample = 500,
  n_samples = 8,
  sigma = 1.2,
  sample_effect = 0,
  lib_loc = 7.5,
  seed = 12345
)

hdw_sim <- sim$seurat
hdw_sim_modules <- merge(
  sim$modules,
  sim$design,
  by = "module",
  sort = FALSE
)
hdw_sim_modules <- hdw_sim_modules[order(hdw_sim_modules$module, -hdw_sim_modules$loading), ]
rownames(hdw_sim_modules) <- NULL

save(hdw_sim, file = "data/hdw_sim.rda", compress = "xz")
save(hdw_sim_modules, file = "data/hdw_sim_modules.rda", compress = "xz")

message("hdw_sim: ", nrow(hdw_sim), " genes x ", ncol(hdw_sim), " cells, ",
        round(file.size("data/hdw_sim.rda") / 1e6, 2), " MB")
