
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
#' @param gene_mapping a dataframe containing gene name mappings between the query and the reference dataset. One column should have the gene name in the query dataset, and anotehr column should have the corresponding gene name in the reference dataset.
#' @param genome_ref_col the column name containing the gene names for the reference dataset
#' @param genome_query_col the column name containing the gene names for the query dataset
#' @param return_raw if TRUE, returns the module preservation statistics, else returns seurat_obj with the stats added to the hdWGCNA experiment.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param wgcna_name_ref The name of the hdWGCNA experiment in the seurat_ref@misc slot
#' @details
#' ModulePreservation performs a statistical test to assess the preservation of co-expression modules identified 
#' in one dataset in an independent dataset. This method is originally described by Langfelder et al. 
#' in the 2011 paper "Is My Network Module Preserved and Reproducible?". This method can be used to 
#' assess biological differences in networks, as well as technical differences / reproducibility across different 
#' batches. 
#'
#' @import Seurat
#' @export
ModulePreservation <- function(
  seurat_obj,
  seurat_ref,
  name,
  n_permutations = 500,
  parallel = FALSE,
  seed = 12345,
  gene_mapping = NULL,
  genome_ref_col = NULL,
  genome_query_col = NULL,
  return_raw = FALSE,
  wgcna_name = NULL,
  wgcna_name_ref = NULL,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  if(is.null(wgcna_name_ref)){wgcna_name_ref <- seurat_ref@misc$active_wgcna}
  CheckWGCNAName(seurat_ref, wgcna_name_ref)
  
  # get datExpr for reference and query:
  datExpr_ref <- GetDatExpr(seurat_ref, wgcna_name_ref)
  datExpr_query <- GetDatExpr(seurat_obj, wgcna_name)

  # change the gene names to match:
  if(!is.null(gene_mapping)){
    gene_match <- match(colnames(datExpr_query), gene_mapping[,genome_query_col])
    gene_mapping <- na.omit(gene_mapping[gene_match,])
    colnames(datExpr_query)  <- gene_mapping[,genome_ref_col]
  }

  genes_keep <- intersect(colnames(datExpr_ref), colnames(datExpr_query))
  datExpr_ref <- datExpr_ref[,genes_keep]
  datExpr_query <- datExpr_query[,genes_keep]

  # set up multiExpr:
  setLabels <- c("ref", "query")
  multiExpr <- list(
    ref = list(data=datExpr_ref),
    query = list(data=datExpr_query)
  )
  ref_modules <- GetModules(seurat_ref, wgcna_name=wgcna_name_ref)
  ref_modules <- ref_modules[genes_keep,]
  ref_modules <- list(ref = ref_modules$module)

  print('ref:')
  print(dim(multiExpr$ref$data))
  print('query:')
  print(dim(multiExpr$query$data))

  # print('Run Module Preservation')
  print(length(ref_modules))

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


#' ModulePreservationNetRep
#'
#' Computes module preservation statistics in a query dataset for a given reference dataset
#'
#' @param seurat_query A Seurat object
#' @param seurat_ref A Seurat object serving as the reference for the module preservation analysis
#' @param name The name to give the module preservation analysis.
#' @param n_permutations Number of permutations for the module preservation test.
#' @param n_threads Number of parallel threads to use for the module preservation test
#' @param TOM_use  The name of the hdWGCNA experiment containing the TOM that will be used for plotting
#' @param gene_mapping a dataframe containing gene name mappings between the query and the reference dataset. One column should have the gene name in the query dataset, and anotehr column should have the corresponding gene name in the reference dataset.
#' @param genome_ref_col the column name containing the gene names for the reference dataset
#' @param genome_query_col the column name containing the gene names for the query dataset
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param wgcna_name_ref The name of the hdWGCNA experiment in the seurat_ref@misc slot
#' @details
#' ModulePreservationNetRep performs a statistical test to assess the preservation of co-expression modules identified 
#' in one dataset in an independent dataset using the NetRep implementation rather than the WGCNA implementation. 
#' This method is originally described by Langfelder et al. in the 2011 paper 
#' "Is My Network Module Preserved and Reproducible?". The NetRep implementation of the module preservation test 
#' is described by Ritchie et al. in the 2016 paper "A Scalable Permutation Approach Reveals Replication and Preservation Patterns of Networks Modules in Large Datasets"
#' This method can be used to assess biological differences in networks, as well as technical differences / reproducibility across different 
#' batches. Note that the outputs of ModulePreservation and ModulePreservationNetRep are not identical. This function requires separate 
#' installation of the NetRep package, which is not directly included as a dependency in hdWGCNA.
#' 
#' @return A Seurat Object containing the module preservation statistics from NetRep.
#' @import Seurat
#' @export
ModulePreservationNetRep <- function(
  seurat_query,
  seurat_ref,
  name,
  n_permutations = 10000,
  n_threads=8,
  gene_mapping = NULL,
  genome_ref_col = NULL,
  genome_query_col = NULL,
  TOM_use = NULL,
  wgcna_name = NULL,
  wgcna_name_ref = NULL,
  ...
){

  # This package is not required for the rest of hdWGCNA so we don't make it a formal dependency
  if (!require("NetRep")) {
    print('Missing package: NetRep')
    print('Installing package: NetRep')
    install.packages('NetRep')
  }

  if(is.null(wgcna_name)){wgcna_name <- seurat_query@misc$active_wgcna}
  CheckWGCNAName(seurat_query, wgcna_name)
  if(is.null(TOM_use)){TOM_use <- wgcna_name}

  if(is.null(wgcna_name_ref)){wgcna_name_ref <- seurat_ref@misc$active_wgcna}
  CheckWGCNAName(seurat_ref, wgcna_name_ref)
    
  # get the modules from the reference dataset
  modules <- GetModules(seurat_ref, wgcna_name = wgcna_name_ref) %>% subset(module != 'grey')
  module_labels <- as.character(modules$module)
  genes_use_ref <- modules$gene_name
  names(module_labels) <- genes_use_ref

  # get datExpr for reference and query:
  datExpr_ref <- as.matrix(GetDatExpr(seurat_ref, wgcna_name_ref))
  datExpr_query <- as.matrix(GetDatExpr(seurat_query, wgcna_name))

  # get TOMs for reference and query:
  TOM_ref <- GetTOM(seurat_ref, wgcna_name_ref)
  TOM_query <- GetTOM(seurat_query, TOM_use)

  # change the gene names to match:
  # this is totally untested so I should test it
  if(!is.null(gene_mapping)){
    gene_match <- match(colnames(datExpr_query), gene_mapping[,genome_query_col])
    gene_mapping <- na.omit(gene_mapping[gene_match,])
    colnames(datExpr_query)  <- gene_mapping[,genome_ref_col]
  }

  # genes to use in teh query dataset
  genes_use_query <- genes_use_ref[genes_use_ref %in% colnames(TOM_query)]

  datExpr_ref <- datExpr_ref[,genes_use_ref]
  datExpr_query <- datExpr_query[,genes_use_query]
  TOM_ref <- TOM_ref[genes_use_ref, genes_use_ref]
  TOM_query <- TOM_query[genes_use_query, genes_use_query]

  # calculate adjacency:
  cor_mat_ref <- WGCNA::adjacency(datExpr_ref, type = "unsigned", power = 1)
  cor_mat_query <- WGCNA::adjacency(datExpr_query, type = "unsigned", power = 1)

  # set up data lists:
  n_list <- list(TOM_ref, TOM_query)
  d_list <- list(datExpr_ref, datExpr_query)
  c_list <- list(cor_mat_ref, cor_mat_query)
  names(n_list) <- c(wgcna_name_ref, wgcna_name)
  names(d_list) <- c(wgcna_name_ref, wgcna_name)
  names(c_list) <- c(wgcna_name_ref, wgcna_name)

  # Assess the preservation of modules in the test dataset.
  preservation_test <- NetRep::modulePreservation(
      network=n_list, 
      data=d_list, 
      correlation=c_list, 
      moduleAssignments=module_labels, 
      discovery=wgcna_name_ref, 
      test=wgcna_name, 
      nPerm=n_permutations, 
      nThreads=n_threads
  )

  # add the result to the seurat object
  seurat_query <- SetModulePreservation(
    seurat_query, 
    preservation_test, 
    name, 
    wgcna_name
  )

  seurat_query

}
