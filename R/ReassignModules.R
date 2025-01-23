#' ReassignModules
#'
#' Reassigns features to different modules
#'
#' @return seurat_obj with the an updated modules table for the selected wgcna experiment
#'
#' @param seurat_obj A Seurat object
#' @param harmonized logical indicating whether or not to use harmonized MEs
#' @param features character vector containing features for manual module reassignment,
#' @param new_modules character vector containing modules to reassign the genes
#' @param ignore logical indicating whether or not to ignore error message about reassigning non-grey features
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' ReassignModules reassigns features with negative kMEs in their assigned module to the
#' module that had the highest kME for that feature. Alternatively, this function
#' can manually assign features to different modules, which can be helpful if
#' certain genes of interest are assigned to the grey module. We generally do not
#' advise reassigning features to new modules if they are in non-grey modules,
#' but the user can do so at their own risk by setting ignore=TRUE.
#'
#' @import Seurat
#' @export
ReassignModules <- function(
  seurat_obj,
  harmonized=TRUE,
  features=NULL,
  new_modules=NULL,
  ignore=FALSE,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get modules and MEs from seurat object
  modules <- GetModules(seurat_obj, wgcna_name)
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  mod_colors <- dplyr::select(modules, c(module, color)) %>% dplyr::distinct()
  mod_cp <- mod_colors$color; names(mod_cp) <- as.character(mod_colors$module)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)

  if(!is.null(features)){
 
    # check validity of input features
    if(!all(features %in% genes_use)){
      stop('Some features are not found in GetWGCNAGenes(seurat_obj).')
    }

     ##############################
    # Manual reassignment
    ##############################

    if(!is.null(new_modules)){

      # check validity of modules
      if(!all(new_modules %in% levels(modules$module))){
        stop('Some module names in new_modules are invalid. Features must be reassigned to existing modules found in GetModules(seurat_obj)')
      }

      # get original modules
      orig_mods <- subset(modules, gene_name %in% features) %>% .$module %>% as.character()
      if(!all(orig_mods == 'grey') & !ignore){
        stop('Attempting to reassign non-grey genes to new modules. Proceed with caution. If you wish to reassign these genes, re-run this function and set ignore=TRUE')
      }

      # set up table for new module assignments
      reassign_df <- data.frame(
        gene_name = as.character(features),
        module = as.character(new_modules)
      )
      reassign_df$module <- factor(as.character(reassign_df$module), levels = levels(modules$module))

      # new colors:
      reassign_df$color <- as.character(mod_cp[as.character(reassign_df$module)])

      # reassign modules and colors
      modules[reassign_df$gene_name,'module'] <- reassign_df$module
      modules[reassign_df$gene_name,'color'] <- reassign_df$color

      # set the modules table
      seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)
      return(seurat_obj)

    }
  } else{

    ##############################
    # reassignment by kME
    ##############################

    # get genes with negative kME in their assigned module:
    neg_df <- do.call(rbind, lapply(mods, function(cur_mod){
      cur <- subset(modules, module == cur_mod)
      cur <- cur[,c('gene_name', 'module', paste0('kME_', cur_mod))]
      names(cur)[3] <- 'kME'
      cur %>% subset(kME < 0)
    }))
    if(nrow(neg_df) == 0){
      return(seurat_obj)
    }
    rownames(neg_df) <- 1:nrow(neg_df)
    features <- neg_df$gene_name

  }

  # get just the kME values from the modules table
  kMEs <- modules[,4:ncol(modules)]
  kMEs <- kMEs[,colnames(kMEs) != "kME_grey"]

  # for each gene with negative kME values in the assigned module,
  # identify the module that had the highest kME
  reassigned <- sapply(features, function(cur_gene){
    cur_kMEs <- kMEs[cur_gene,]
    max_kME <- max(cur_kMEs)
    if(max_kME < 0){
      return('kME_grey')
    }
    colnames(kMEs)[which(cur_kMEs == max_kME)]
  })

  # add the reassigned modules to the modules table
  reassigned <- do.call(rbind, strsplit(reassigned, 'kME_'))[,2]
  reassigned <- factor(as.character(reassigned), levels=levels(modules$module))

  # new colors:
  reassigned_colors <- as.character(mod_cp[as.character(reassigned)])

  # reassign modules and colors
  modules[features,'module'] <- reassigned
  modules[features,'color'] <- reassigned_colors



  # set the modules table
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)
  seurat_obj
}



#' ResetModuleNames
#'
#' Reset the uname of each hdWGCNA module
#'
#' @param seurat_obj A Seurat object
#' @param new_name string containing the base name to re-name the modules
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
ResetModuleNames <- function(
  seurat_obj,
  new_name = "M",
  reset_levels = FALSE,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  old_mods <- levels(modules$module)
  if('grey' %in% modules$module){
    nmods <- length(old_mods) - 1
  } else{
      nmods <- length(old_mods)
  }

  # if it's a named list:
  if(class(new_name) == 'list'){
    if(all(names(new_name) %in% old_mods)){
      ix <- match(names(new_name), old_mods)
      new_names <- old_mods
      new_names[ix] <- as.character(new_name)
      new_names <- new_names[new_names != 'grey']
    } else{
      stop("Some module names present in new_name are not found in this hdWGCNA experiment.")
    }
  } else if(length(new_name) == 1){
    new_names <- paste0(new_name, 1:nmods)
  } else if(length(new_name) == nmods){
    new_names <- new_name
  } else{
    stop("Invalid input for new_name.")
  }

  if('grey' %in% modules$module){

    grey_ind <- which(old_mods == 'grey')

    # account for when grey is first / last
    if(grey_ind == 1){
      new_names <- c('grey', new_names)
    } else if(grey_ind == length(old_mods)){
      new_names <- c(new_names, 'grey')
    } else{
      new_names <- c(new_names[1:(grey_ind-1)], 'grey', new_names[grey_ind:length(new_names)])
    }

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
    me_colnames <- colnames(hMEs)
    ix <- match(me_colnames, new_mod_df$old)
    colnames(hMEs) <- new_mod_df$new[ix]
    seurat_obj <- SetMEs(seurat_obj, hMEs, harmonized=TRUE, wgcna_name)
  }

  # update ME table
  MEs <- GetMEs(seurat_obj, harmonized=FALSE, wgcna_name)
  if(!is.null(MEs)){
    me_colnames <- colnames(MEs)
    ix <- match(me_colnames, new_mod_df$old)
    colnames(MEs) <- new_mod_df$new[ix]
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

  # update degree 
  degree_df <- GetDegrees(seurat_obj, wgcna_name)
  if(!is.null(degree_df)){
    degree_df$module <- factor(
      new_mod_df[match(degree_df$module, new_mod_df$old),'new'],
      levels = as.character(new_mod_df$new)
    )
    seurat_obj <- SetDegrees(seurat_obj, degree_df, wgcna_name)
  }
  
  # return the seurat object
  seurat_obj

}


############################
# Reset module names:
###########################

#' ResetModuleColors
#'
#' Reset the unique color for each hdWGCNA module
#'
#' @param seurat_obj A Seurat object
#' @param new_colors a character vector containing the new colors
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
ResetModuleColors <- function(
  seurat_obj,
  new_colors,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  mod_colors_df <- dplyr::select(modules, c(module, color)) %>%
    dplyr::distinct() %>% dplyr::arrange(module)
  mod_colors <- mod_colors_df$color
  if('grey' %in% modules$mod){
    grey_ind <- which(mod_colors == 'grey')
  } else{
    grey_ind <- NA
  }

  # case where we give a named list:
  if(class(new_colors) == 'list'){
    ix <- match(names(new_colors), mod_colors_df$module)
    mod_colors_df[ix, 'color'] <- as.character(new_colors)
    new_colors <- mod_colors_df$color
  } else{
    if(is.na(grey_ind)){
      new_colors <- new_colors
    } else if(grey_ind == 1){
      new_colors <- c('grey', new_colors)
    } else if(grey_ind == length(mod_colors)){
      new_colors <- c(new_colors, 'grey')
    } else{
      new_colors <- c(new_colors[1:(grey_ind-1)], 'grey', new_colors[grey_ind:length(new_colors)])
    }
  }

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

  # update hdWGCNA dendrogram:
  net <- GetNetworkData(seurat_obj, wgcna_name)
  net$colors <- new_color_df[match(net$colors, new_color_df$old),'new']
  seurat_obj <- SetNetworkData(seurat_obj, net, wgcna_name)

  seurat_obj

}
