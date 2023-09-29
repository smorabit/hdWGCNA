#' ComputeModuleEigengene
#'
#' Internal helper function that computes module eigengene for a single module.
#'
#' @param seurat_obj A Seurat object
#' @param cur_mod name of a module found in seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_net$colors
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are "linear", "poisson", or "negbinom"
#' @param pc_dim Which PC to use as the module eigengene? Default to 1.
#' @param assay Assay in seurat_obj to compute module eigengenes. Default is DefaultAssay(seurat_obj)
#' @param wgcna_name name of the WGCNA experiment
#' @import Seurat
#' @import harmony
ComputeModuleEigengene <- function(
  seurat_obj,
  cur_mod,
  modules,
  group.by.vars=NULL,
  verbose=TRUE,
  vars.to.regress = NULL,
  scale.model.use = 'linear',
  pc_dim = 1,
  assay = NULL,
  wgcna_name=NULL, ...
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get the assay
  if(is.null(assay)){assay <- DefaultAssay(seurat_obj)}
  if(dim(seurat_obj@assays[[assay]]@data)[1] == 0){
    stop(paste0("Normalized data slot not found in selected assay ", assay))
  }

  # get genes in this module:
  cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name %>% as.character()

  # subset seurat object by these genes only:
  X_dat <- GetAssayData(seurat_obj, slot='data', assay = assay)[cur_genes,]
  if(dim(seurat_obj@assays[[assay]]@counts)[1] == 0){
    X <- X_dat
  } else{
    X <- GetAssayData(seurat_obj, slot='counts', assay = assay)[cur_genes,]
  }

  # create seurat obj with just these genes
  cur_seurat <- CreateSeuratObject(X, assay = assay, meta.data = seurat_obj@meta.data)
  cur_seurat <- SetAssayData(cur_seurat, slot='data', new.data=X_dat, assay=assay)

  # scale the subsetted expression dataset:
  if(is.null(vars.to.regress)){
    cur_seurat <- ScaleData(cur_seurat, features=rownames(cur_seurat), model.use=scale.model.use)
  } else if(all(vars.to.regress %in% colnames(seurat_obj@meta.data))){
    cur_seurat <- ScaleData(cur_seurat, features=rownames(cur_seurat), model.use=scale.model.use, vars.to.regress=vars.to.regress)
  } else{
    stop(paste0("Some variables specified in vars.to.regress are not found in seurat_obj@meta.data"))
  }

  # compute average expression of each gene
  cur_expr <- GetAssayData(cur_seurat, slot='data')
  expr <- Matrix::t(cur_expr)
  averExpr <- Matrix::rowSums(expr) / ncol(expr)

  # run PCA with Seurat function
  cur_pca <- Seurat::RunPCA(
    cur_seurat,
    features = cur_genes,
    reduction.key=paste0('pca', cur_mod),
    verbose=verbose, ...
  )@reductions$pca
  pc <- cur_pca@cell.embeddings[,pc_dim]
  pc_loadings <- cur_pca@feature.loadings[,pc_dim]

  # correlate average expression with eigengene
  pca_cor <- cor(averExpr, pc)

  # run harmony
  if(!is.null(group.by.vars)){

    # add this PCA as its own reduction in the seurat object
    seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
      embeddings = cur_pca@cell.embeddings,
      assay = assay
    )

    cur_harmony <- harmony::RunHarmony(
      seurat_obj,
      group.by.vars=group.by.vars,
      reduction.use="ME", verbose=verbose, assay.use=assay, ...
    )@reductions$harmony
    ha <- cur_harmony@cell.embeddings[,pc_dim]
    ha_loadings <- cur_pca@feature.loadings[,pc_dim]

    if(pca_cor < 0){
      cur_harmony@cell.embeddings[,pc_dim] <- -ha
      ha_loadings <- -ha_loadings
    }

    # add harmonized PCA as its own reduction in the seurat object
    seurat_obj@reductions$ME_harmony <- Seurat::CreateDimReducObject(
      embeddings = cur_harmony@cell.embeddings,
      assay = assay
    )

    seurat_obj <- SetMELoadings(
      seurat_obj,
      loadings=ha_loadings,
      harmonized=TRUE,
      wgcna_name=wgcna_name
    )

  }

  if(pca_cor < 0){
    cur_pca@cell.embeddings[,pc_dim] <- -pc
    pc_loadings <- -pc_loadings
  }

  # add this PCA as its own reduction in the seurat object
  seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
    embeddings = cur_pca@cell.embeddings,
    assay = assay
  )

  seurat_obj <- SetMELoadings(
    seurat_obj,
    loadings=pc_loadings,
    harmonized=FALSE,
    wgcna_name=wgcna_name
  )

  # return seurat object
  seurat_obj

}

#' ModuleEigengenes
#'
#' Computes module eigengenes for co-expression modules
#'
#' @return seurat_obj with the module eigengenes computed for the selected wgcna experiment
#'
#' @param seurat_obj A Seurat object
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are "linear", "poisson", or "negbinom"
#' @param pc_dim Which PC to use as the module eigengene? Default to 1.
#' @param assay Assay in seurat_obj to compute module eigengenes. Default is DefaultAssay(seurat_obj)
#' @param exclude_grey logical determining whether to compute MEs for the grey module
#' @param wgcna_name name of the WGCNA experiment
#'
#' @details
#' ModuleEigengenes summarizes the gene expression signatures of entire co-expression
#' modules. This is done by performing singular value decomposition (SVD) on a
#' subset of the scaled expression matrix containing only features assigned to each module.
#' The module eigengene (ME), defined as the first dimension of the SVD matrix, retains the most variation, and we use
#' this vector as a summary of gene expression for the whole module.
#'
#' The module gene expression matrix is first scaled using the Seurat ScaleData
#' function. The user can optionally adjust for covariates of interest in this step using the
#' vars.to.regress parameter. Additionally, the module eigengenes themselves can
#' be adjusted for technical biases such as sequencing batch, dataset of origin,
#' or other factors using the Harmony algorithm with the group.by.vars parameter.
#'
#' @import Seurat
#' @export
ModuleEigengenes <- function(
  seurat_obj,
  group.by.vars=NULL,
  modules=NULL,
  vars.to.regress = NULL,
  scale.model.use = 'linear',
  verbose=TRUE,
  assay = NULL,
  pc_dim = 1,
  exclude_grey = FALSE,
  wgcna_name=NULL, ...
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  if(!CheckWGCNAName(seurat_obj, wgcna_name)){
    stop(paste0("Invalid wgcna_name supplied: ", wgcna_name))
  }  

  # are we going to run Harmony?
  harmonized = !is.null(group.by.vars)

  if(harmonized & !any(grepl("ScaleData", seurat_obj@commands))){
      stop('Need to run ScaleData before running ModuleEigengenes with group.by.vars option.')
  }

  # get the assay
  if(is.null(assay)){assay <- DefaultAssay(seurat_obj)}

  # check for the data slot in this assay
  if(dim(seurat_obj@assays[[assay]]@data)[1] == 0){
    stop(paste0("Normalized data slot not found in selected assay ", assay))
  }

  # exclude grey doesn't work yet:
  if(exclude_grey){exclude_grey <- FALSE}


  me_list <- list()
  harmonized_me_list <- list()

  # re-set feature loadings:
  seurat_obj <- SetMELoadings(
    seurat_obj,
    loadings=c(""),
    harmonized=FALSE,
    wgcna_name=wgcna_name
  )
  if(harmonized){
    seurat_obj <- SetMELoadings(
      seurat_obj,
      loadings=c(""),
      harmonized=TRUE,
      wgcna_name=wgcna_name
    )
  }

  # get modules from Seurat object, else use provided modules
  if(is.null(modules)){
    modules <- GetModules(seurat_obj, wgcna_name)
    projected <- FALSE
  } else{
    projected <- TRUE
  }

  # get list of modules:
  mods <- levels(modules$module)
  mods_loop <- mods

  # exclude grey:
  if(exclude_grey){
    mods_loop <- mods[mods != 'grey']
  }

  # loop over modules:
  for(cur_mod in mods_loop){

    print(cur_mod)

    # compute module eigengenes for this module
    seurat_obj <- ComputeModuleEigengene(
      seurat_obj = seurat_obj,
      cur_mod = cur_mod,
      modules=modules,
      group.by.vars=group.by.vars,
      vars.to.regress = vars.to.regress,
      scale.model.use = scale.model.use,
      verbose=verbose,
      pc_dim = pc_dim,
      assay = assay,
      wgcna_name=wgcna_name,
      ...
    )

    # add module eigengene to ongoing list
    cur_me <- seurat_obj@reductions$ME@cell.embeddings[,pc_dim]
    me_list[[cur_mod]] <- cur_me

    # run harmony
    if(harmonized){
      # add module eigengene to ongoing list
      cur_harmonized_me <- seurat_obj@reductions$ME_harmony@cell.embeddings[,pc_dim]
      harmonized_me_list[[cur_mod]] <- cur_harmonized_me
    }
  }

  # merge module eigengene lists into a dataframe, order modules, add to Seurat obj
  me_df <- do.call(cbind, me_list)
  if(!projected){me_df <- WGCNA::orderMEs(me_df)}
  seurat_obj <- SetMEs(seurat_obj, me_df, harmonized=FALSE, wgcna_name)

  # merge harmonized module eigengene lists into a dataframe, add to Seurat obj
  if(!is.null(group.by.vars)){
    hme_df <- do.call(cbind, harmonized_me_list)
    if(!projected){hme_df <- WGCNA::orderMEs(hme_df)}
    seurat_obj <- SetMEs(seurat_obj, hme_df, harmonized=TRUE, wgcna_name)
  }

  # set module factor levels based on order
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  modules$module <- factor(
    as.character(modules$module),
    levels=mods
  )
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

  # remove temp dim reductions by setting to NULL
  seurat_obj@reductions$ME <- NULL
  seurat_obj@reductions$ME_harmony <- NULL

  # return seurat object
  seurat_obj
}

#' ModuleExprScore
#'
#' Computes a module score for each co-expression module using Seurat AddModuleScore or UCell.
#'
#' @return seurat_obj with module scores computed for the selected wgcna experiment
#'
#' @param seurat_obj A Seurat object
#' @param n_genes the number of genes to use for each module, ranked by kME. Setting n_genes = 'all' uses all of the genes in a module
#' @param method selected method for module scoring, valid choices are "Seurat" or "UCell"
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' ModuleExprScore provides an alternative function to ModuleEigengenes for summarizing
#' the expression level of each module. The user can choose between Seurat AddModuleScore
#' or UCell using the method parameter.
#'
#' @export
ModuleExprScore <- function(
  seurat_obj,
  n_genes = 25,
  method='Seurat',
  wgcna_name=NULL,
  ...
 ){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  if(!CheckWGCNAName(seurat_obj, wgcna_name)){
    stop(paste0("Invalid wgcna_name supplied: ", wgcna_name))
  }  

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # use all genes in the module?
  if(n_genes == "all"){
    gene_list <- lapply(mods, function(cur_mod){
      subset(modules, module == cur_mod) %>% .$gene_name
    })
  } else{
    gene_list <- lapply(mods, function(cur_mod){
      cur <- subset(modules, module == cur_mod)
      cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
        top_n(n_genes) %>% .$gene_name
    })
  }
  names(gene_list) <- mods

  # compute Module Scores with Seurat or UCell:
  if(method == "Seurat"){
    mod_scores <- Seurat::AddModuleScore(
      seurat_obj, features=gene_list, ...
    )@meta.data
  } else if(method == "UCell"){
    mod_scores <- UCell::AddModuleScore_UCell(
      seurat_obj, features=gene_list, ...
    )@meta.data
  } else{
    stop("Invalid method selection. Valid choices are Seurat, UCell")
  }

  mod_scores <- mod_scores[,(ncol(mod_scores)-length(mods)+1):ncol(mod_scores)]

  # rename module scores:
  colnames(mod_scores) <- mods

  # order cols
  col_order <- levels(modules$module)
  col_order <- col_order[col_order != 'grey']
  mod_scores <- mod_scores[,col_order]

  # add module scores to seurat object
  seurat_obj <- SetModuleScores(seurat_obj, mod_scores, wgcna_name)

  seurat_obj

}


#' AverageModuleExpr
#'
#' Computes module eigengenes for all WGCNA co-expression modules
#'
#' @param seurat_obj A Seurat object
#' @export
AvgModuleExpr <- function(seurat_obj, n_genes = 25, wgcna_name=NULL, ...){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  if(!CheckWGCNAName(seurat_obj, wgcna_name)){
    stop(paste0("Invalid wgcna_name supplied: ", wgcna_name))
  }  
  
  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get wgcna params
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  # datExpr for full expression dataset
  datExpr <- GetAssayData(
    seurat_obj,
    assay=params$metacell_assay,
    slot=params$metacell_slot
  )

  # use all genes in the module?
  if(n_genes == "all"){
    gene_list <- lapply(mods, function(cur_mod){
      subset(modules, module == cur_mod) %>% .$gene_name
    })
  } else{
    gene_list <- lapply(mods, function(cur_mod){
      cur <- subset(modules, module == cur_mod)
      cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
        top_n(n_genes) %>% .$gene_name
    })
  }
  names(gene_list) <- mods

  # for each module, compute average expression of genes:
  avg_mods <- t(do.call(rbind,lapply(gene_list, function(cur_genes){colMeans(datExpr[cur_genes,])})))

  # order cols
  col_order <- levels(modules$module)
  col_order <- col_order[col_order != 'grey']
  avg_mods <- avg_mods[,col_order]

  # add avg module expression to Seurat object
  seurat_obj <- SetAvgModuleExpr(seurat_obj, avg_mods, wgcna_name)
  seurat_obj
}
