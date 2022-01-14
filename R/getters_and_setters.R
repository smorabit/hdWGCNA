
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


#' SetDatExpr
#'
#' This function sets up the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' SetDatExpr(pbmc)
SetDatExpr <- function(seurat_obj, use_metacells=TRUE, wgcna_name=NULL, group.by=NULL, group_name=NULL, ...){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get parameters from seurat object
  params <- GetWGCNAParams(seurat_obj, wgcna_name)
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  assay <- params$metacell_assay

  # use metacells or whole seurat object?
  if(use_metacells){
    s_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  } else{
    s_obj <- seurat_obj
  }

  # columns to group by
  if(!is.null(group.by)){
    cells <- s_obj@meta.data %>% subset(get(group.by) == group_name) %>% rownames
  } else{
    cells <- colnames(s_obj)
  }

  # get expression data from seurat obj
  datExpr <- as.data.frame(
    Seurat::GetAssayData(
      s_obj,
      assay=assay,
      slot='data'
    )[genes_use,cells]
  )

  # transpose data
  datExpr <- as.data.frame(t(datExpr))

  # only get good genes:
  gene_list = GetWGCNAGenes(seurat_obj)[WGCNA::goodGenes(datExpr, ...)]

  # update the WGCNA gene list:
  seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)

  datExpr <- datExpr[,gene_list]

  # set the datExpr in the Seurat object
  seurat_obj@misc[[wgcna_name]]$datExpr <- datExpr

  # return seurat obj
  seurat_obj
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
  hMEs <- GetMEs(seurat_obj, wgcna_name)
  colnames(hMEs) <- new_mod_df$new
  seurat_obj <- SetMEs(seurat_obj, hMEs, harmonized=TRUE, wgcna_name)

  # update ME table
  MEs <- GetMEs(seurat_obj, harmonized=FALSE, wgcna_name)
  colnames(MEs) <- new_mod_df$new
  seurat_obj <- SetMEs(seurat_obj, MEs, harmonized=FALSE, wgcna_name)

  # update module scores:
  module_scores <- GetModuleScores(seurat_obj, wgcna_name)
  if(!("grey" %in% colnames(module_scores))){
    colnames(module_scores) <- new_mod_df$new[new_mod_df$new != 'grey']
  } else {
    colnames(module_scores) <- new_mod_df$new
  }
  seurat_obj <- SetModuleScores(seurat_obj, module_scores, wgcna_name)

  # update average module expression:
  avg_exp <- GetAvgModuleExpr(seurat_obj, wgcna_name)
  if(!("grey" %in% colnames(avg_exp))){
    colnames(avg_exp) <- new_mod_df$new[new_mod_df$new != 'grey']
  } else {
    colnames(avg_exp) <- new_mod_df$new
  }
  seurat_obj <- SetAvgModuleExpr(seurat_obj, avg_exp, wgcna_name)

  # update enrichr table:
  enrich_table <- GetEnrichrTable(seurat_obj, wgcna_name)
  enrich_table$module <- factor(
    new_mod_df[match(enrich_table$module, new_mod_df$old),'new'],
    levels = as.character(new_mod_df$new)
  )
  seurat_obj <- SetEnrichrTable(seurat_obj, enrich_table, wgcna_name)

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
  mod_colors <- select(modules, c(module, color)) %>%
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
    old = mod_colors ,
    new = new_colors
  )

  modules$color <- new_color_df[match(modules$color, new_color_df$old),'new']

  # set module table
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

  seurat_obj

}
