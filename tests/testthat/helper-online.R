# Opt-in tests for functions that call web services or need very large
# annotation packages. They never run on CRAN or in the regular CI jobs;
# run them by setting the environment variable HDWGCNA_ONLINE_TESTS=true and
# calling devtools::test(filter = 'online'), or trigger the 'online-tests' GitHub workflow by hand.
skip_if_no_online_tests <- function() {
  testthat::skip_on_cran()
  if (!identical(tolower(Sys.getenv('HDWGCNA_ONLINE_TESTS')), 'true')) {
    testthat::skip('Online tests are off; set HDWGCNA_ONLINE_TESTS=true to run them')
  }
  testthat::skip_if_offline()
}

# Gene sets with unambiguous Enrichr hits, used to give simulated modules
# real gene symbols.
online_gene_sets <- list(
  ribosome = c(
    'RPL3', 'RPL4', 'RPL5', 'RPL6', 'RPL7', 'RPL8', 'RPL10', 'RPL11', 'RPL13',
    'RPL14', 'RPL15', 'RPL18', 'RPL19', 'RPL21', 'RPL23', 'RPL27', 'RPS2',
    'RPS3', 'RPS4X', 'RPS5', 'RPS6', 'RPS8', 'RPS9', 'RPS11', 'RPS12', 'RPS14',
    'RPS15', 'RPS16', 'RPS18', 'RPS19'
  ),
  cell_cycle = c(
    'CDK1', 'CCNB1', 'CCNB2', 'CCNA2', 'CDC20', 'BUB1', 'BUB1B', 'PLK1',
    'AURKA', 'AURKB', 'MAD2L1', 'CDC25C', 'CCNE1', 'CDK2', 'MCM2', 'MCM3',
    'MCM4', 'MCM5', 'MCM6', 'MCM7', 'PCNA', 'E2F1', 'TOP2A', 'MKI67',
    'KIF11', 'CENPF', 'TTK', 'ESPL1', 'CDKN3', 'CHEK1'
  )
)

# A finished simulated network whose first two modules carry the gene
# symbols above, so Enrichr returns meaningful terms for them.
sim_wgcna_real_genes <- function() {
  so <- sim_wgcna('kmes')
  modules <- GetModules(so)
  mods <- setdiff(levels(modules$module), 'grey')[1:2]
  for (i in 1:2) {
    ix <- which(modules$module == mods[i])
    n <- min(length(ix), length(online_gene_sets[[i]]))
    ix <- ix[order(-modules[ix, paste0('kME_', mods[i])])][seq_len(n)]
    modules$gene_name <- as.character(modules$gene_name)
    modules$gene_name[ix] <- online_gene_sets[[i]][seq_len(n)]
  }
  SetModules(so, modules)
}
