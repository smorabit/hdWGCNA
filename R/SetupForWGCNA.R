

#' SetupForWGCNA
#'
#' Create a slot in a Seurat object to store hdWGCNA data
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name name of the WGCNA experiment
#' @param features list of features to use for WGCNA
#' @param metacell_location name of the WGCNA experiment to copy the metacell object from
#' @param ... additional parameters to pass to SelectNetworkGenes
#' 
#' @details 
#' SetupForWGCNA creates a new slot in the Seurat object (seurat_obj) to store an hdWGCNA experiment
#' with a given name (wgcna_name). This function calls on SelectNetworkGenes to specify which features
#' will be used for network analysis. If there is another hdWGCNA experiment already in the seurat_obj,
#' the same metacell/metaspot object can be used by specifying the name of the hdWGCNA experiment 
#' to the metacell_location parameter.
#'
#' @export
SetupForWGCNA <- function(
  seurat_obj, wgcna_name,
  features = NULL,
  metacell_location = NULL,
  ...
){

  # set the active WGCNA variable
  seurat_obj <- SetActiveWGCNA(seurat_obj, wgcna_name)

  # select genes for WGCNA, else use pre-selected:
  if(is.null(features)){
    seurat_obj <- SelectNetworkGenes(seurat_obj, wgcna_name=wgcna_name, ...)
  } else{

    # set the selected genes:
   # seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list=features, wgcna_name=wgcna_name)
    seurat_obj <- SelectNetworkGenes(
        seurat_obj, 
        wgcna_name=wgcna_name, 
        gene_select = 'custom', 
        gene_list = features
      )

  }

  # give the user a warning if there's too few genes (under 200?)
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name=wgcna_name)
  if(length(genes_use) < 500){
    warning(paste0(length(genes_use), ' features selected. You may wish to change the settings to include more features in your analysis.'))
  }

  if(!is.null(metacell_location)){
    if(!(metacell_location %in% names(seurat_obj@misc))){
      stop('metacell_location not found in seurat_obj@misc')
    }
    seurat_obj <- SetMetacellObject(seurat_obj, metacell_location, wgcna_name)
  }

  seurat_obj
}
