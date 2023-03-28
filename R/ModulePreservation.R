
#' ModulePreservation
#'
#' Computes module preservation statistics in a query dataset for a given reference dataset
#'
#' @param seurat_obj A Seurat object
#' @param seurat_ref A Seurat object serving as the reference for the module preservation analysis
#' @param name The name to give the module preservation analysis.
#' @param n_permutations Number of permutations for the module preservation test.
#' @param parallel logical determining whether to run preservation analysis in parallel
#' @param seed random seed for the permutation analysis.
#' @param return_raw if TRUE, returns the module preservation statistics, else returns seurat_obj with the stats added to the hdWGCNA experiment.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param wgcna_name_ref The name of the hdWGCNA experiment in the seurat_ref@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModulePreservation
ModulePreservation <- function(
  seurat_obj,
  seurat_ref,
  name,
  n_permutations = 500,
  parallel = FALSE,
  seed = 12345,
  gene_mapping = NULL,
  genome1_col = NULL,
  genome2_col = NULL,
  return_raw = FALSE,
  wgcna_name = NULL,
  wgcna_name_ref = NULL,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  if(is.null(wgcna_name_ref)){wgcna_name_ref <- seurat_ref@misc$active_wgcna}

  print('Setup datasets')

  # get datExpr for reference and query:
  datExpr_ref <- GetDatExpr(seurat_ref, wgcna_name_ref)
  datExpr_query <- GetDatExpr(seurat_obj, wgcna_name)

  # change the gene names to match:
  if(!is.null(gene_mapping)){
    gene_match <- match(colnames(datExpr_query), gene_mapping[,genome2_col])
    gene_mapping <- na.omit(gene_mapping[gene_match,])
    colnames(datExpr_query)  <- gene_mapping[,genome1_col]

    print(head(colnames(datExpr_query)))
    print(head(GetModules(seurat_obj)$gene_name))
  }

  # set up multiExpr:
  setLabels <- c("ref", "query")
  multiExpr <- list(
    ref = list(data=datExpr_ref),
    query = list(data=datExpr_query)
  )
  ref_modules <- list(ref = GetModules(seurat_ref)$module)

  print('ref:')
  print(dim(datExpr_ref))
  print('query:')
  print(dim(datExpr_query))

  print('Run Module Preservation')
  print(ref_modules)

  # run the module preservation test:
  mp <- WGCNA::modulePreservation(
    multiExpr,
    ref_modules,
    referenceNetworks = 1,
    nPermutations = n_permutations,
    randomSeed = seed,
    quickCor = 0,
    parallelCalculation = parallel,
    ...
  )

  if(return_raw){return(mp)}

  # get the stats obs and stats Z tables
  ref <- 1; test <- 2;
  statsObs <- cbind(
    mp$quality$observed[[ref]][[test]],
    mp$preservation$observed[[ref]][[test]]
  )
  statsZ <- cbind(
    mp$quality$Z[[ref]][[test]],
    mp$preservation$Z[[ref]][[test]]
  )

  # add stats to the seurat object:
  mod_pres <- list('obs'  = statsObs, 'Z' = statsZ)
  seurat_obj <- SetModulePreservation(seurat_obj, mod_pres, name, wgcna_name)

  seurat_obj

}
