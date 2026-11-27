############################
# TF Match Matrix (not stored within the WGCNA slot)
###########################

#' SetMotifMatrix
#'
#' @param seurat_obj A Seurat object
#' @param tf_match matrix containing tf-promoter matches
#' @keywords scRNA-seq
#' @export
SetMotifMatrix <- function(seurat_obj, tf_match) {
  # make a spot for the motif info if it's not already there:
  if (is.null(seurat_obj@misc$motifs)) {
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$tf_match_matrix <- tf_match
  seurat_obj
}


#' GetMotifMatrix
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
GetMotifMatrix <- function(seurat_obj) {
  seurat_obj@misc$motifs$tf_match_matrix
}


############################
# Motif table
###########################

#' SetMotifs
#'
#' @param seurat_obj A Seurat object
#' @param motif_df dataframe containing info about the motifs being analyzed
#' @keywords scRNA-seq
#' @export
SetMotifs <- function(seurat_obj, motif_df) {
  # make a spot for the motif info if it's not already there:
  if (is.null(seurat_obj@misc$motifs)) {
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$motif_df <- motif_df
  seurat_obj
}


#' GetMotifs
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
GetMotifs <- function(seurat_obj) {
  seurat_obj@misc$motifs$motif_df
}


############################
# PFM List
###########################


#' SetPFMList
#'
#' @param seurat_obj A Seurat object
#' @param pfm_list list of pfm objects
#' @keywords scRNA-seq
#' @export
SetPFMList <- function(seurat_obj, pfm_list) {
  # make a spot for the motif info if it's not already there:
  if (is.null(seurat_obj@misc$motifs)) {
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$pfm_list <- pfm_list
  seurat_obj
}


#' GetPFMList
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
GetPFMList <- function(seurat_obj) {
  seurat_obj@misc$motifs$pfm_list
}


############################
# TF Target Genes:
###########################

#' SetMotifTargets
#'
#' @param seurat_obj A Seurat object
#' @param motif_targets list of motifs and their target genes
#' @keywords scRNA-seq
#' @export
SetMotifTargets <- function(seurat_obj, motif_targets) {
  # make a spot for the motif info if it's not already there:
  if (is.null(seurat_obj@misc$motifs)) {
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$motif_targets <- motif_targets
  seurat_obj
}


#' GetMotifTargets
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
GetMotifTargets <- function(seurat_obj) {
  seurat_obj@misc$motifs$motif_targets
}


############################
# motif overlap
###########################


#' SetMotifOverlap
#'
#' @param seurat_obj A Seurat object
#' @param overlap_df dataframe containing motif-module overlap info
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMotifOverlap <- function(seurat_obj, overlap_df, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$motif_module_overlaps <- overlap_df
  seurat_obj
}


#' GetMotifOverlap
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMotifOverlap <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$motif_module_overlaps
}


############################
# motif scores
###########################


#' SetMotifScores
#'
#' @param seurat_obj A Seurat object
#' @param tf_scores dataframe of tf motif target scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMotifScores <- function(seurat_obj, tf_scores, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$motif_target_scores <- tf_scores
  seurat_obj
}


#' GetMotifScores
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMotifScores <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$motif_target_scores
}


############################
# getters and setters
###########################

#' SetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param tf_net dataframe storing the TF network info in ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFNetwork <- function(seurat_obj, tf_net, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_net <- tf_net
  seurat_obj
}


#' GetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFNetwork <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_net
}


#' SetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param tf_eval dataframe storing the TF network evaluation info from ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFEval <- function(seurat_obj, tf_eval, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_eval <- tf_eval
  seurat_obj
}


#' GetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFEval <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_eval
}


#' SetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param tf_regulons dataframe storing the TF regulon info from AssignTFRegulons
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFRegulons <- function(seurat_obj, tf_regulons, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_regulons <- tf_regulons
  seurat_obj
}


#' GetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFRegulons <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_regulons
}




#' SetRegulonScores
#'
#' @param seurat_obj A Seurat object
#' @param regulon_scores dataframe storing the TF regulon scores
#' @param target_type dataframe storing the TF regulon scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetRegulonScores <- function(seurat_obj, regulon_scores, target_type, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)


  # if regulon scores have not been set, make a list to store them
  if (is.null(seurat_obj@misc[[wgcna_name]]$regulon_scores)) {
    tmp <- list(regulon_scores)
    names(tmp) <- target_type
    seurat_obj@misc[[wgcna_name]]$regulon_scores <- tmp
  } else {
    seurat_obj@misc[[wgcna_name]]$regulon_scores[[target_type]] <- regulon_scores
  }
  seurat_obj
}


#' GetRegulonScores
#'
#' @param seurat_obj A Seurat object
#' @param target_type dataframe storing the TF regulon scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetRegulonScores <- function(seurat_obj, target_type, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get the regulon scores
  seurat_obj@misc[[wgcna_name]]$regulon_scores[[target_type]]
}


#' SetEnrichRegulonTable
#'
#' @param seurat_obj A Seurat object
#' @param enrich_table dataframe storing the results of running enrichr
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetEnrichrRegulonTable <- function(seurat_obj, enrich_table, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set enrichr table
  seurat_obj@misc[[wgcna_name]]$enrichr_regulon_table <- enrich_table
  seurat_obj
}


#' GetEnrichrRegulonTable
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetEnrichrRegulonTable <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$enrichr_regulon_table
}
