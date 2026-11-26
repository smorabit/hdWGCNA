# Runtime / memory benchmarks of the core hdWGCNA pipeline on simulated data.
#
# Usage (from the package root; one scale per process so peak RSS is per scale):
#   /usr/bin/time -l pixi run Rscript bench/run_benchmarks.R <scale> [label]
#   e.g.  /usr/bin/time -l pixi run Rscript bench/run_benchmarks.R S baseline
# or run all scales with bench/run_all.sh
#
# Output: bench/results/<label>_<scale>.csv with one row per pipeline step:
#   step, seconds (wall), r_mem_max_mb (peak R heap during the step, from gc()),
#   plus dataset size and environment info.
#
# Notes
# - r_mem_max_mb only covers memory allocated by R; memory allocated in C/C++/Fortran
#   (e.g. inside WGCNA) is invisible to gc(). Process-level peak RSS comes from
#   /usr/bin/time -l (see run_all.sh, which appends it to the CSV).
# - Threads: WGCNA threads are fixed via allowWGCNAThreads(n_threads).

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})
source("data-raw/sim_functions.R")

args <- commandArgs(trailingOnly = TRUE)
scale <- if (length(args) >= 1) args[1] else "S"
label <- if (length(args) >= 2) args[2] else format(Sys.Date())
n_threads <- as.integer(Sys.getenv("HDW_BENCH_THREADS", "4"))

scales <- list(
  S = list(n_genes = 1000, cells_per_sample = 500,  n_samples = 8),   # = shipped hdw_sim
  M = list(n_genes = 2000, cells_per_sample = 2000, n_samples = 8),   # 16k cells
  L = list(n_genes = 4000, cells_per_sample = 2000, n_samples = 16)   # 32k cells
)
stopifnot(scale %in% names(scales))
cfg <- scales[[scale]]

WGCNA::allowWGCNAThreads(n_threads)

results <- list()
timed <- function(step, expr) {
  invisible(gc(reset = TRUE))
  t <- system.time(value <- force(expr))[["elapsed"]]
  g <- gc()
  mem <- sum(g[, ncol(g)])  # "max used" (Mb) for Ncells + Vcells
  results[[length(results) + 1]] <<- data.frame(step = step, seconds = t, r_mem_max_mb = mem)
  message(sprintf("%-22s %8.1f s  %8.0f MB", step, t, mem))
  value
}

sim <- timed("simulate", do.call(simulate_hdwgcna_data, c(cfg, list(sigma = 1.2, lib_loc = 7.5))))
so <- sim$seurat
rm(sim)
so <- timed("preprocess", preprocess_sim(so))

so <- timed("SetupForWGCNA",
  SetupForWGCNA(so, gene_select = "fraction", fraction = 0.05, wgcna_name = "bench"))
so <- timed("MetacellsByGroups",
  MetacellsByGroups(so, group.by = c("cell_type", "Sample"), reduction = "pca",
                    k = 20, max_shared = 10, min_cells = 40, ident.group = "cell_type",
                    verbose = FALSE))
so <- timed("NormalizeMetacells", NormalizeMetacells(so))
so <- timed("SetDatExpr",
  SetDatExpr(so, group_name = "A", group.by = "cell_type", assay = "RNA", layer = "data"))
so <- timed("TestSoftPowers", suppressWarnings(TestSoftPowers(so, networkType = "signed")))
tom_dir <- tempfile("TOM")
so <- timed("ConstructNetwork",
  ConstructNetwork(so, tom_outdir = tom_dir, minModuleSize = 20, verbose = 0))
so <- timed("ModuleEigengenes", ModuleEigengenes(so, group.by.vars = "Sample", verbose = FALSE))
so <- timed("ModuleConnectivity",
  ModuleConnectivity(so, group.by = "cell_type", group_name = "A"))

out <- do.call(rbind, results)
out$scale <- scale
out$n_cells <- ncol(so)
out$n_genes <- nrow(so)
out$n_network_genes <- length(GetWGCNAGenes(so))
out$n_metacells_A <- nrow(GetDatExpr(so))
out$threads <- n_threads
out$blas <- basename(extSoftVersion()[["BLAS"]])
out$r_version <- paste(R.version$major, R.version$minor, sep = ".")
out$hdwgcna_commit <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                               error = function(e) NA)
out$date <- format(Sys.time(), "%Y-%m-%d %H:%M")

dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
f <- sprintf("bench/results/%s_%s.csv", label, scale)
utils::write.csv(out, f, row.names = FALSE)
message("written: ", f)
