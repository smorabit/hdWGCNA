#' ModuleConnectivity
#'
#' Computes eigengene-based connectivity (kME) based on module eigengenes and the
#' feature expression in selected cell populations / spatial regions.
#'
#' @return seurat_obj with the kMEs computed for the selected wgcna experiment
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing grouping info, ie clusters or celltypes
#' @param group_name name of the group(s) in group.by to use for kME calculation
#' @param corFnc character string specifying the function to be used to calculate co-expression similarity. Defaults to bicor. Any function returning values
#' @param corOptions 	character string specifying additional arguments to be passed to the function given by corFnc. Use "use = 'p', method = 'spearman'" to obtain Spearman correlation. Use "use = 'p'" to obtain Pearson correlation.
#' @param harmonized logical determining whether to use harmonized MEs for kME calculation
#' @param assay Assay in seurat_obj containing expression information.
#' @param slot Slot in seurat_obj, default to normalized 'data' slot.
#' @param sparse logical indicating whether or not to run the correlation using a sparse matrix.
#' @param reassign_modules logical indicating whether or not to reassign genes to different co-expression modules if their kME value in the assigned module is negative.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' ModuleConnectivity computes the eigengene-based connectivity (kME) of each feature
#' and each module. This is done by correlating the gene expression signal with
#' the module eigengenes for each module. Features with high kME values have greater
#' network connectivity within that module, and we can identify module hub genes
#' by taking the top genes for each module by kME values.
#'
#' We recommend using the group.by and group_name parameters that were previously
#' used in the SetDatExpr function, so only the relevant cell populations / spatial
#' regions are considered for computing kMEs.
#'
#' Depending on the parameters for network analysis, sometimes a feature has a negative kME
#' in the module that it was originally assigned. This discrepancy arises since
#' the network construction is done on the metacell/metaspot matrix but the kME
#' calculation is done on the cell/spot matrix. By default, we account for this
#' by reassigning features with negative kMEs in their assigned module to the
#' module that had the highest kME for that feature.
#'
#' @import WGCNA
#' @import Seurat
#' @import Matrix
#' @export
ModuleConnectivity <- function(
  seurat_obj,
  group.by = NULL,
  group_name = NULL,
  corFnc = 'bicor',
  corOptions = "use='p'",
  harmonized=TRUE,
  assay = NULL,
  slot = 'data',
  layer = 'data',
  sparse = TRUE,
  reassign_modules = TRUE,
  wgcna_name = NULL,
  ...
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  
  # get module df, wgcna genes, and wgcna params:
  modules <- GetModules(seurat_obj, wgcna_name)
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  genes_use <- as.character(modules$gene_name)
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  if(is.null(assay)){assay <- DefaultAssay(seurat_obj)}

  if(!is.null(group.by)){
    cells.use <- seurat_obj@meta.data %>% subset(get(group.by) %in% group_name) %>% rownames
    MEs <- MEs[cells.use,]
  } else{
    cells.use <- colnames(seurat_obj)
  }

  # get expression data:
  if(CheckSeurat5()){
    exp_mat <- SeuratObject::LayerData(
      seurat_obj, assay=assay, layer=layer
    )[genes_use,cells.use]
  } else{
    exp_mat <- GetAssayData(
      seurat_obj,
      assay=assay,
      slot=slot
    )[genes_use,cells.use]
  }

  if(sparse){
    
    kMEs <- corSparse(
      X = Matrix::t(exp_mat),
      Y = as.matrix(MEs)
    )
    rownames(kMEs) <- genes_use
    kMEs <- as.data.frame(kMEs)
  } else{
    datExpr <- t(as.matrix(exp_mat))
    kMEs <- WGCNA::signedKME(
      datExpr,
      MEs,
      outputColumnName = "kME",
      corFnc = corFnc,
      corOptions = corOptions,
      ...
    )

  }

  # add module color to the kMEs table
  modules <- modules[,1:3]
  mods <- levels(modules$module)
  colnames(kMEs) <- colnames(MEs)
  kMEs <- kMEs[,mods]
  colnames(kMEs) <- paste0("kME_", colnames(kMEs))
  kMEs <- cbind(modules, kMEs)

  # update the modules table in the Seurat object:
  seurat_obj <- SetModules(seurat_obj, kMEs, wgcna_name)

  # reassign modules for genes with negative kME values
  if(reassign_modules & length(mods) > 2){
    seurat_obj <- ReassignModules(seurat_obj, harmonized=harmonized, wgcna_name=wgcna_name)
  }

  seurat_obj

}
