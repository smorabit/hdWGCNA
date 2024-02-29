
#' FindAllDMEs
#'
#' Function to compare expression levels of co-expression modules between two sets of cell barcodes.
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing cell grouping information
#' @param features indicate whether to use "MEs" or "ModuleScores" for the comparison
#' @param harmonized logical determining whether or not to use harmonized MEs
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param add_missing logical determining whether or not to add missing modules back into the resulting dataframe with NA values.
#' @param ... Additional parameters for the FindMarkers function
#' @keywords scRNA-seq
#' @export
#' @return A dataframe contaning differential ME results
#' @examples
#' FindAllDMEs
FindAllDMEs <- function(
  seurat_obj,
  group.by,
  features = 'MEs',
  harmonized=TRUE,
  add_missing=FALSE,
  wgcna_name=NULL,
  test.use='wilcox',
  only.pos=FALSE,
  logfc.threshold = 0,
  min.pct=0,
  verbose=FALSE,
  pseudocount.use=0,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get list of groups
  groups <- as.character(unique(seurat_obj@meta.data[[group.by]]))

  # what features are we using?
  if(features == 'MEs'){
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  } else if(features == 'ModuleScores'){
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else{
    stop('Invalid selection for features. Valid choices are MEs or ModuleScores.')
  }

  # remove grey module
  MEs <- MEs[, colnames(MEs) != 'grey']

  # set all negative values to zero
  MEs[MEs < 0] <- 0

  # transpose MEs so cells are columns and rows are MEs
  MEs <- t(MEs)

  # create new assay for MEs:
  ME_assay <- Seurat::CreateAssayObject(MEs)

  DMEs_list <- list()
  for(cur_group in groups){

    print(cur_group)

    barcodes1 <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] == cur_group]
    barcodes2 <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] != cur_group]

    # run differential test on the ME assay
    DMEs <- FindMarkers(
      ME_assay,
      cells.1 = barcodes1,
      cells.2 = barcodes2,
      slot='counts',
      test.use=test.use,
      only.pos=only.pos,
      logfc.threshold=logfc.threshold,
      min.pct=min.pct,
      verbose=verbose,
      pseudocount.use=pseudocount.use,
      ...
    )
    DMEs$module <- rownames(DMEs)
    DMEs$group <- cur_group

    # add missing modules if specified
    if(add_missing){
      missing_mods <- rownames(MEs)[!(rownames(MEs) %in% DMEs$module)]
      for(cur_mod in missing_mods){
        DMEs[cur_mod,] <- NA
        DMEs[cur_mod,'module'] <- cur_mod
      }
    }

    # add to ongoing list:
    DMEs_list[[cur_group]] <- DMEs

  }

  DMEs <- do.call(rbind, DMEs_list)
  DMEs
}


#' FindDMEs
#'
#' Function to compare expression levels of co-expression modules between two sets of cell barcodes.
#'
#' @param seurat_obj A Seurat object
#' @param barcodes1 character vector containing cell barcodes for the first group to test. Positive fold-change means up-regulated in this group.
#' @param barcodes2 character vector containing cell barcodes for the second group to test. Negative fold-change means up-regulated in this group.
#' @param features indicate whether to use "MEs" or "ModuleScores" for the comparison
#' @param harmonized logical determining whether or not to use harmonized MEs
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param add_missing logical determining whether or not to add missing modules back into the resulting dataframe with NA values.
#' @param ... Additional parameters for the FindMarkers function
#' @keywords scRNA-seq
#' @export
#' @return A dataframe contaning differential ME results
#' @examples
#' FindDMEs
FindDMEs <- function(
  seurat_obj,
  barcodes1,
  barcodes2,
  features = 'MEs',
  harmonized=TRUE,
  wgcna_name=NULL,
  add_missing=FALSE,
  test.use='wilcox',
  only.pos=FALSE,
  logfc.threshold = 0,
  min.pct=0,
  verbose=FALSE,
  pseudocount.use=0,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # ensure that selected barcodes are in the seurat obj
  if(!(all(barcodes1 %in% colnames(seurat_obj)))){
    stop('Some barcodes in barcodes1 not found in colnames(seurat_obj).')
  }
  if(!(all(barcodes2 %in% colnames(seurat_obj)))){
    stop('Some barcodes in barcodes2 not found in colnames(seurat_obj).')
  }

  # check for overlap in the two groups of bacrodes:
  if(length(intersect(barcodes1, barcodes2)) > 0){
    stop('Some barcodes overlap in barcodes1 and barcodes2')
  }

  # what features are we using?
  if(features == 'MEs'){
    MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  } else if(features == 'ModuleScores'){
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else{
    stop('Invalid selection for features. Valid choices are MEs or ModuleScores.')
  }

  # remove grey module
  MEs <- MEs[, colnames(MEs) != 'grey']

  # set all negative values to zero
  MEs[MEs < 0] <- 0

  # transpose MEs so cells are columns and rows are MEs
  MEs <- t(MEs)

  # create new assay for MEs:
  ME_assay <- Seurat::CreateAssayObject(MEs)

  # run differential test on the ME assay
  DMEs <- FindMarkers(
    ME_assay,
    cells.1 = barcodes1,
    cells.2 = barcodes2,
    slot='counts',
    test.use=test.use,
    only.pos=only.pos,
    logfc.threshold=logfc.threshold,
    min.pct=min.pct,
    verbose=verbose,
    pseudocount.use=pseudocount.use,
    ...
  )
  DMEs$module <- rownames(DMEs)

  # add missing modules if specified
  if(add_missing){
    missing_mods <- rownames(MEs)[!(rownames(MEs) %in% DMEs$module)]
    for(cur_mod in missing_mods){
      DMEs[cur_mod,] <- NA
      DMEs[cur_mod,'module'] <- cur_mod
    }
  }

  DMEs

}
