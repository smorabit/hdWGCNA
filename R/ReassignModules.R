#' ReassignModules
#'
#' Reassigns features to different modules
#'
#' @return seurat_obj with the an updated modules table for the selected wgcna experiment
#'
#' @param seurat_obj A Seurat object
#' @param harmonized logical indicating whether or not to use harmonized MEs
#' @param features character vector containing features for manual module reassignment,
#' @param new_modules character vector containing modules to ,
#' @param ignore logical indicating whether or not to ignore error message about reassigning non-grey features
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' ReassignModules reassigs features with negative kMEs in their assigned module to the
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

  # get modules and MEs from seurat object
  modules <- GetModules(seurat_obj, wgcna_name)
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct()
  mod_cp <- mod_colors$color; names(mod_cp) <- as.character(mod_colors$module)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)


  if(!is.null(features) & !is.null(new_modules)){

    ##############################
    # Manual reassignment
    ##############################

    # check validity of input features
    if(!all(features %in% genes_use)){
      stop('Some features are not found in GetWGCNAGenes(seurat_obj).')
    }

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
      stop('No genes to reassign, all kMEs of assigned modules are greater than 0.')
    }
    rownames(neg_df) <- 1:nrow(neg_df)


    # get just the kME values from the modules table
    kMEs <- modules[,4:ncol(modules)]
    kMEs <- kMEs[,colnames(kMEs) != "kME_grey"]

    # for each gene with negative kME values in the assigned module,
    # identify the module that had the highest kME
    reassigned <- sapply(neg_df$gene_name, function(cur_gene){
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
    modules[neg_df$gene_name,'module'] <- reassigned
    modules[neg_df$gene_name,'module'] <- reassigned_colors

  }

  # set the modules table
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)
  seurat_obj
}
