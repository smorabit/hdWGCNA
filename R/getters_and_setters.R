
############################
# Active WGCNA
###########################

SetActiveWGCNA <- function(seurat_obj, wgcna_name){

  # set the active_wgcna variable
  seurat_obj@misc$active_wgcna <- wgcna_name

  # initialize empty list for this WGCNA if it doesn't exist yet
  if(!(wgcna_name %in% names(seurat_obj@misc))){
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]] <- list()
  }
  seurat_obj
}

GetActiveWGCNA <- function(seurat_obj){
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]
}

# get any WGCNA data, but by default get the active
GetWGCNA <- function(seurat_obj, wgcna_name=NULL){

  # test if wgcna_name is valid (TODO)

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  seurat_obj@misc[[wgcna_name]]
}

############################
# WGCNA Group
###########################

SetWGCNAGroup <- function(seurat_obj, group, wgcna_name){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # add gene list to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_group <- group
  seurat_obj
}

GetWGCNAGroup <- function(seurat_obj, wgcna_name){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_group
}


############################
# metacell object
###########################

SetMetacellObject <- function(seurat_obj, metacell_obj, wgcna_name=NULL){
  # test if wgcna_name is valid (TODO)
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # add metacell obj to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj <- metacell_obj
  seurat_obj
}

GetMetacellObject <- function(seurat_obj,  wgcna_name=NULL){
  # test if wgcna_name is valid (TODO)

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj
}

############################
# WGCNA genes
###########################

SetWGCNAGenes <- function(seurat_obj, gene_list, wgcna_name=NULL){

  # test if wgcna_name is valid (TODO)
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # add gene list to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_genes <- gene_list
  seurat_obj
}

GetWGCNAGenes <- function(seurat_obj, wgcna_name=NULL){
  # test if wgcna_name is valid (TODO)
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_genes
}

############################
# datExpr
###########################


#' spr
#'
#' This function sets up the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @param group_name A string containing a group present in the provided group.by column or in the Seurat Idents.
#' @param use_metacells A logical determining if we use the metacells (TRUE) or the full expression matrix (FALSE)
#' @param group.by A string containing the name of a column in the Seurat object with cell groups (clusters, cell types, etc). If NULL (default), scWGCNA uses the Seurat Idents as the group.
#' @param multi.group.by A string containing the name of a column in the Seurat object with groups for consensus WGCNA (dataset, sample, condition, etc)
#' @param multi_group_name A string containing the name of a group present in the multi.group.by column.
#' @param wgcna_name A string containing the name of the WGCNA slot in seurat_obj@misc. Default = NULL which retrieves the currently active WGCNA data
#' @keywords scRNA-seq
#' @export
#' @examples
#' SetDatExpr(pbmc)
SetDatExpr <- function(
  seurat_obj,
  group_name,
  use_metacells=TRUE,
  group.by=NULL,
  multi.group.by = NULL,
  multi_group_name = NULL,
  return_seurat = TRUE,
  wgcna_name=NULL,
  slot = 'data',
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get parameters from seurat object
  params <- GetWGCNAParams(seurat_obj, wgcna_name)
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  assay <- params$metacell_assay

  print('n_genes:')
  print(length(genes_use))
  print(head(genes_use))

  # use metacells or whole seurat object?
  if(use_metacells){
    s_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  } else{
    s_obj <- seurat_obj
  }

  # get the metadata from the seurat object:
  seurat_meta <- s_obj@meta.data
  print(dim(seurat_meta))

  # columns to group by for cluster/celltype
  if(!is.null(group.by)){
    seurat_meta <- seurat_meta %>% subset(get(group.by) == group_name)
  }
  print(dim(seurat_meta))

  # subset further if multiExpr:
  if(!is.null(multi.group.by)){
    seurat_meta <- seurat_meta %>% subset(get(multi.group.by) == multi_group_name)
  }
  print(dim(seurat_meta))

  # get list of cells to use
  cells <- rownames(seurat_meta)
  print('cells:')
  print(head(cells))
  print(length(cells))

  # get expression data from seurat obj
  datExpr <- as.data.frame(
    Seurat::GetAssayData(
      s_obj,
      assay=assay,
      slot=slot
    )[genes_use,cells]
  )

  # transpose data
  datExpr <- as.data.frame(t(datExpr))

  print(dim(datExpr))

  # only get good genes:
  if(is.null(multi.group.by)){
    gene_list = genes_use[WGCNA::goodGenes(datExpr, ...)]
    datExpr <- datExpr[,gene_list]
  }

  print(dim(datExpr))

  if(return_seurat){

    # update the WGCNA gene list:
    seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)

    # set the datExpr in the Seurat object
    seurat_obj@misc[[wgcna_name]]$datExpr <- datExpr
    out <- seurat_obj
  } else{
    out <- datExpr
  }
  out
}


#' GetDatExpr
#'
#' This function gets the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' GetDatExpr(pbmc)
GetDatExpr <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$datExpr

}


#' SetMultiExpr
#'
#' This function sets up the expression matrix input for consensus WGCNA based on
#' the metacell expression matrix, or from the full expression matrix.
#'
#' @param seurat_obj A Seurat object
#' @param group_name A string containing a group present in the provided group.by column or in the Seurat Idents.
#' @param use_metacells A logical determining if we use the metacells (TRUE) or the full expression matrix (FALSE)
#' @param group.by A string containing the name of a column in the Seurat object with cell groups (clusters, cell types, etc). If NULL (default), scWGCNA uses the Seurat Idents as the group.
#' @param multi.group.by A string containing the name of a column in the Seurat object with groups for consensus WGCNA (dataset, sample, condition, etc)
#' @param multi_groups A character vecrtor containing the names of
#' @param wgcna_name A string containing the name of the WGCNA slot in seurat_obj@misc. Default = NULL which retrieves the currently active WGCNA data
#' @keywords scRNA-seq
#' @export
#' @examples
#' SetDatExpr(pbmc)
SetMultiExpr <- function(
  seurat_obj,
  group_name,
  use_metacells=TRUE,
  group.by=NULL,
  multi.group.by = NULL,
  multi_groups = NULL,
  wgcna_name=NULL,
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get the WGCNA genes:
  gene_names <- GetWGCNAGenes(seurat_obj, wgcna_name)

  # get the different groups present if not specified by the user:
  if(is.null(multi_groups)){
    multi_groups <- unique(seurat_obj@meta.data[[multi.group.by]])
  } else{
    seurat_groups <- unique(seurat_obj@meta.data[[multi.group.by]])
    if(sum(multi_groups %in% seurat_groups) != length(multi_groups)){
      stop('Some or all groups specified in multi_groups not found in seurat_obj@meta.data[,multi.group.by]')
    }
  }

  # get the datExpr for each group
  datExpr_list <- lapply(multi_groups, function(x){
    SetDatExpr(
      seurat_obj,
      group_name = group_name,
      group.by = group.by,
      multi.group.by = multi.group.by,
      multi_group_name = x,
      return_seurat = FALSE,
      wgcna_name = wgcna_name
    ) %>% as.matrix
  })
  datExpr_list

  # convert to multiExpr, get good genes:
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  genes_use <- WGCNA::goodGenesMS(multiExpr)
  gene_names <- gene_names[genes_use]

  # subset the multiExpr by the good genes::
  datExpr_list <- lapply(1:length(multiExpr), function(i){
    multiExpr[[i]]$data[,genes_use]
  })
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  names(multiExpr) <- multi_groups

  # update the WGCNA gene list:
  seurat_obj <- SetWGCNAGenes(seurat_obj, gene_names, wgcna_name)

  # set the multiExpr in the Seurat object
  seurat_obj@misc[[wgcna_name]]$multiExpr <- multiExpr
  seurat_obj

}


#' GetMultiExpr
#'
#' This function gets the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' GetDatExpr(pbmc)
GetMultiExpr <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$multiExpr

}

############################
# WGCNA params
###########################

SetWGCNAParams <- function(seurat_obj, params, wgcna_name=NULL){

  if(is.null(GetActiveWGCNA(seurat_obj)$wgcna_params)){
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params <- params
  } else{
    for(i in 1:length(params)){
      seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[[names(params)[i]]] <- params[[i]]
    }
  }
  seurat_obj
}

GetWGCNAParams <- function(seurat_obj, wgcna_name=NULL){
  # test if wgcna_name is valid (TODO)

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_params
}

############################
# SoftPower Table
###########################

SetPowerTable <- function(seurat_obj, power_table, wgcna_name=NULL){
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # add power table to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_powerTable <- power_table
  seurat_obj
}

GetPowerTable <- function(seurat_obj, wgcna_name=NULL){
  # get data from active assay if wgcna_name is not given

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_powerTable
}

############################
# WGCNA Network
###########################

SetNetworkData <- function(seurat_obj, net, wgcna_name=NULL){
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # add network data to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_net <- net
  seurat_obj
}

GetNetworkData <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_net
}

############################
# WGCNA modules dataframe
###########################

SetModules <- function(seurat_obj, mod_df, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # set module df
  seurat_obj@misc[[wgcna_name]]$wgcna_modules <- mod_df
  seurat_obj
}


GetModules <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_modules
}


############################
# Module Eigengenes
###########################

SetMEs <- function(seurat_obj, MEs, harmonized=TRUE, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # harmonized MEs?
  if(harmonized){
    seurat_obj@misc[[wgcna_name]]$hMEs <- MEs
  } else{
    seurat_obj@misc[[wgcna_name]]$MEs <- MEs
  }
  seurat_obj
}

#' GetMEs
#'
#' Function to retrieve module eigengens from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#'  ModuleEigengenes(pbmc)
GetMEs <- function(seurat_obj, harmonized=TRUE, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get harmonized MEs?
  if(harmonized == TRUE && !is.null(seurat_obj@misc[[wgcna_name]]$hMEs)){
    MEs <- seurat_obj@misc[[wgcna_name]]$hMEs
  } else{
    MEs <- seurat_obj@misc[[wgcna_name]]$MEs
  }
  MEs
}

############################
# GO term table
###########################

SetEnrichrTable <- function(seurat_obj, enrich_table, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # set enrichr table
  seurat_obj@misc[[wgcna_name]]$enrichr_table <- enrich_table
  seurat_obj
}


GetEnrichrTable <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$enrichr_table
}


############################
# Module Scores
###########################

SetModuleScores <- function(seurat_obj, mod_scores, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$module_scores <- mod_scores
  seurat_obj
}


GetModuleScores <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$module_scores
}

############################
# Average Module Expression
###########################

SetAvgModuleExpr <- function(seurat_obj, avg_mods, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$avg_modules <- avg_mods
  seurat_obj
}


GetAvgModuleExpr <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$avg_modules
}

############################
# TOM
###########################

GetTOM <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get WGCNA genes:
  gene_names <- GetWGCNAGenes(seurat_obj, wgcna_name)

  # load TOM
  tom_files <- GetNetworkData(seurat_obj, wgcna_name)$TOMFiles
  load(tom_files[[1]])

  TOM <- as.matrix(consTomDS)
  rownames(TOM) <- gene_names; colnames(TOM) <- gene_names
  TOM

}

############################
# ROC Curve
###########################

SetROCData <- function(seurat_obj, roc_info, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$roc_data <- roc_info
  seurat_obj
}


GetROCData <- function(seurat_obj,  wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$roc_data
}


############################
# TF Match Matrix (not stored within the WGCNA slot)
###########################

SetMotifMatrix <- function(seurat_obj, tf_match){

  # make a spot for the motif info if it's not already there:
  if(is.null(seurat_obj@misc$motifs)){
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$tf_match_matrix <- tf_match
  seurat_obj
}

GetMotifMatrix <- function(seurat_obj){
  seurat_obj@misc$motifs$tf_match_matrix
}


############################
# Motif table
###########################

SetMotifs <- function(seurat_obj, motif_df){

  # make a spot for the motif info if it's not already there:
  if(is.null(seurat_obj@misc$motifs)){
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$motif_df <- motif_df
  seurat_obj
}

GetMotifs <- function(seurat_obj){
  seurat_obj@misc$motifs$motif_df
}

############################
# PFM List
###########################

SetPFMList <- function(seurat_obj, pfm_list){

  # make a spot for the motif info if it's not already there:
  if(is.null(seurat_obj@misc$motifs)){
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$pfm_list <- pfm_list
  seurat_obj
}

GetPFMList <- function(seurat_obj){
  seurat_obj@misc$motifs$pfm_list
}


############################
# TF Target Genes:
###########################

SetMotifTargets <- function(seurat_obj, motif_targets){

  # make a spot for the motif info if it's not already there:
  if(is.null(seurat_obj@misc$motifs)){
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$motif_targets <- motif_targets
  seurat_obj
}

GetMotifTargets <- function(seurat_obj){
  seurat_obj@misc$motifs$motif_targets
}


############################
# motif overlap
###########################

SetMotifOverlap <- function(seurat_obj, overlap_df, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$motif_module_overlaps <- overlap_df
  seurat_obj
}

GetMotifOverlap <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$motif_module_overlaps
}


############################
# motif scores
###########################

SetMotifScores <- function(seurat_obj, tf_scores, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$motif_target_scores <- tf_scores
  seurat_obj
}

GetMotifScores <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$motif_target_scores
}

############################
# ModuleUMAP
###########################

SetModuleUMAP <- function(seurat_obj, umap_df, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$module_umap <- umap_df
  seurat_obj
}

GetModuleUMAP <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$module_umap
}

############################
# ModuleTraitCorrelation
###########################

SetModuleTraitCorrelation <- function(seurat_obj, mt_cor, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$mt_cor <- mt_cor
  seurat_obj
}

GetModuleTraitCorrelation <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$mt_cor
}


############################
# Reset module names:
###########################

ResetModuleNames <- function(
  seurat_obj,
  new_name = "M",
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  old_mods <- levels(modules$module)

  new_names <- paste0(new_name, 1:(length(old_mods)-1))
  grey_ind <- which(old_mods == 'grey')

  # account for when grey is first / last
  if(grey_ind == 1){
    new_names <- c('grey', new_names)
  } else if(grey_ind == length(old_mods)){
    new_names <- c(new_names, 'grey')
  } else{
    new_names <- c(new_names[1:(grey_ind-1)], 'grey', new_names[grey_ind:length(new_names)])
  }

  # update kMEs
  new_kMEs <- paste0('kME_', new_names)
  colnames(modules) <- c(colnames(modules)[1:3], new_kMEs)

  # update module names
  new_mod_df <- data.frame(
    old = old_mods ,
    new = new_names
  )

  modules$module <- factor(
    new_mod_df[match(modules$module, new_mod_df$old),'new'],
    levels = as.character(new_mod_df$new)
  )

  # set module table
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

  # update hME table:
  hMEs <- GetMEs(seurat_obj, harmonized=TRUE, wgcna_name)
  if(!is.null(hMEs)){
    colnames(hMEs) <- new_mod_df$new
    seurat_obj <- SetMEs(seurat_obj, hMEs, harmonized=TRUE, wgcna_name)
  }

  # update ME table
  MEs <- GetMEs(seurat_obj, harmonized=FALSE, wgcna_name)
  if(!is.null(MEs)){
    colnames(MEs) <- new_mod_df$new
    seurat_obj <- SetMEs(seurat_obj, MEs, harmonized=FALSE, wgcna_name)
  }

  # update module scores:
  module_scores <- GetModuleScores(seurat_obj, wgcna_name)
  if(!is.null(module_scores)){
    if(!("grey" %in% colnames(module_scores))){
      colnames(module_scores) <- new_mod_df$new[new_mod_df$new != 'grey']
    } else {
      colnames(module_scores) <- new_mod_df$new
    }
    seurat_obj <- SetModuleScores(seurat_obj, module_scores, wgcna_name)
  }

  # update average module expression:
  avg_exp <- GetAvgModuleExpr(seurat_obj, wgcna_name)
  if(!is.null(avg_exp)){
    if(!("grey" %in% colnames(avg_exp))){
      colnames(avg_exp) <- new_mod_df$new[new_mod_df$new != 'grey']
    } else {
      colnames(avg_exp) <- new_mod_df$new
    }
    seurat_obj <- SetAvgModuleExpr(seurat_obj, avg_exp, wgcna_name)
  }

  # update enrichr table:
  enrich_table <- GetEnrichrTable(seurat_obj, wgcna_name)
  if(!is.null(enrich_table)){
    enrich_table$module <- factor(
      new_mod_df[match(enrich_table$module, new_mod_df$old),'new'],
      levels = as.character(new_mod_df$new)
    )
    seurat_obj <- SetEnrichrTable(seurat_obj, enrich_table, wgcna_name)
  }

  # update ROC info:
  # THIS DOES NOT UPDATE THE ROC OBJECTS THEMSELVES!!!
  roc_data <- GetROCData(seurat_obj, wgcna_name)
  if(!is.null(roc_data)){
    roc_data$roc$module <- factor(
      new_mod_df[match(roc_data$roc$module, new_mod_df$old),'new'],
      levels = as.character(new_mod_df$new)
    )
    roc_data$conf$module <- factor(
      new_mod_df[match(roc_data$conf$module, new_mod_df$old),'new'],
      levels = as.character(new_mod_df$new)
    )
    seurat_obj <- SetROCData(seurat_obj, roc_data, wgcna_name)
  }

  # update motif overlap
  overlap_df <- GetMotifOverlap(seurat_obj, wgcna_name)
  if(!is.null(overlap_df)){
    overlap_df$module <- factor(
      new_mod_df[match(overlap_df$module, new_mod_df$old),'new'],
      levels = as.character(new_mod_df$new)
    )
    seurat_obj <- SetMotifOverlap(seurat_obj, overlap_df, wgcna_name)
  }

  # update module umap:
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  if(!is.null(umap_df)){
    umap_df$module <- factor(
      new_mod_df[match(umap_df$module, new_mod_df$old),'new'],
      levels = as.character(new_mod_df$new)
    )
    seurat_obj <- SetModuleUMAP(seurat_obj, umap_df, wgcna_name)
  }


  seurat_obj

}


############################
# Reset module names:
###########################

ResetModuleColors <- function(
  seurat_obj,
  new_colors,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  mod_colors <- dplyr::select(modules, c(module, color)) %>%
    distinct %>% arrange(module) %>% .$color
  grey_ind <- which(mod_colors == 'grey')

  if(grey_ind == 1){
    new_colors <- c('grey', new_colors)
  } else if(grey_ind == length(mod_colors)){
    new_colors <- c(new_colors, 'grey')
  } else{
    new_colors <- c(new_colors[1:(grey_ind-1)], 'grey', new_colors[grey_ind:length(new_colors)])
  }
  new_colors

  new_color_df <- data.frame(
    old = mod_colors,
    new = new_colors
  )

  modules$color <- new_color_df[match(modules$color, new_color_df$old),'new']

  # set module table
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)


  # update motif overlap
  overlap_df <- GetMotifOverlap(seurat_obj, wgcna_name)
  if(!is.null(overlap_df)){
    overlap_df$color <- new_color_df[match(overlap_df$color, new_color_df$old),'new']
    seurat_obj <- SetMotifOverlap(seurat_obj, overlap_df, wgcna_name)
  }

  # update module umap:
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  if(!is.null(umap_df)){
    umap_df$color <- new_color_df[match(umap_df$color, new_color_df$old),'new']
    seurat_obj <- SetModuleUMAP(seurat_obj, umap_df, wgcna_name)
  }


  seurat_obj

}
