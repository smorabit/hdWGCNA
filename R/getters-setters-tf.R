############################
# TF Match Matrix (not stored within the WGCNA slot)
###########################

#' SetMotifMatrix
#'
#' @param seurat_obj A Seurat object.
#' @param tf_match A logical gene-by-motif matrix of TF motif matches in gene promoters.
#' @return The input Seurat object with the motif matrix stored in \code{seurat_obj@misc$motifs$tf_match_matrix}.
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
#' @param seurat_obj A Seurat object.
#' @return A logical gene-by-motif matrix of TF motif matches in gene promoters, as computed by \code{MotifScan}.
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
#' @param seurat_obj A Seurat object.
#' @param motif_df A data.frame containing information about the motifs being analyzed.
#' @return The input Seurat object with the motif table stored in \code{seurat_obj@misc$motifs$motif_df}.
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
#' @param seurat_obj A Seurat object.
#' @return A data.frame of motif information (motif ID, motif name, number of targets, and gene name).
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
#' @param seurat_obj A Seurat object.
#' @param pfm_list A list of position frequency matrices (e.g. a PFMatrixList).
#' @return The input Seurat object with the PFM list stored in \code{seurat_obj@misc$motifs$pfm_list}.
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
#' @param seurat_obj A Seurat object.
#' @return The list of position frequency matrices used for the motif scan.
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
#' @param seurat_obj A Seurat object.
#' @param motif_targets A named list of motifs and their target genes.
#' @return The input Seurat object with the motif targets stored in \code{seurat_obj@misc$motifs$motif_targets}.
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
#' @param seurat_obj A Seurat object.
#' @return A named list with one character vector of target genes per motif.
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
#' @param seurat_obj A Seurat object.
#' @param overlap_df A data.frame of motif-module overlap statistics.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the overlap table stored in \code{seurat_obj@misc[[wgcna_name]]$motif_module_overlaps}.
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
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of motif-module overlap statistics, as computed by \code{OverlapModulesMotifs}, or NULL if not computed.
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
#' @param seurat_obj A Seurat object.
#' @param tf_scores A data.frame of TF motif target scores.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the motif target scores stored in \code{seurat_obj@misc[[wgcna_name]]$motif_target_scores}.
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
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of TF motif target scores, or NULL if not set.
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
#' @param seurat_obj A Seurat object.
#' @param tf_net A data.frame of TF-gene importance scores from \code{ConstructTFNetwork}.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the TF network stored in \code{seurat_obj@misc[[wgcna_name]]$tf_net}.
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
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of TF-gene importance scores, as computed by \code{ConstructTFNetwork}.
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
#' @param seurat_obj A Seurat object.
#' @param tf_eval A data.frame of TF network model evaluation results from \code{ConstructTFNetwork}.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the TF network evaluation stored in \code{seurat_obj@misc[[wgcna_name]]$tf_eval}.
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
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of TF network model evaluation results, as computed by \code{ConstructTFNetwork}.
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
#' @param seurat_obj A Seurat object.
#' @param tf_regulons A data.frame of TF regulons from \code{AssignTFRegulons}.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the TF regulons stored in \code{seurat_obj@misc[[wgcna_name]]$tf_regulons}.
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
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of TF regulons, as computed by \code{AssignTFRegulons}.
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
#' @param seurat_obj A Seurat object.
#' @param regulon_scores A data.frame of TF regulon scores (cells by TFs).
#' @param target_type The regulon target type the scores were computed for ("positive", "negative", or "both"); used as the storage key.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the scores stored in \code{seurat_obj@misc[[wgcna_name]]$regulon_scores[[target_type]]}.
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
#' @param seurat_obj A Seurat object.
#' @param target_type The regulon target type to retrieve scores for ("positive", "negative", or "both").
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of TF regulon scores (cells by TFs) for \code{target_type}, as computed by \code{RegulonScores}, or NULL if not computed.
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


#' SetEnrichrRegulonTable
#'
#' @param seurat_obj A Seurat object.
#' @param enrich_table A data.frame of Enrichr results for TF regulons.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the regulon Enrichr table stored in \code{seurat_obj@misc[[wgcna_name]]$enrichr_regulon_table}.
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
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of Enrichr results for TF regulons, as computed by \code{RunEnrichrRegulons}, or NULL if not computed.
#' @keywords scRNA-seq
#' @export
GetEnrichrRegulonTable <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$enrichr_regulon_table
}
