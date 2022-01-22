
#' NormalizeMetacells
#'
#' Wrapper function to run Seurat's NormalizeData function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' NormalizeMetadata
NormalizeMetacells <- function(seurat_obj, wgcna_name=NULL, ...){
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::NormalizeData(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}

#' ScaleMetacells
#'
#' Wrapper function to run Seurat's ScaleData function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' ScaleMetadata
ScaleMetacells <- function(seurat_obj, wgcna_name=NULL, ...){
  if(!exists("features")){
    features = VariableFeatures(seurat_obj)
  }
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::ScaleData(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}

#' RunPCAMetacells
#'
#' Wrapper function to run Seurat's RunPCA function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' NormalizeMetadata
RunPCAMetacells <- function(seurat_obj, wgcna_name=NULL, ...){
  metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  metacell_obj <- Seurat::RunPCA(metacell_obj, ...)
  SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)
}

#' RunHarmonyMetacells
#'
#' Wrapper function to run harmony's RunHarmony function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' NormalizeMetadata
RunHarmonyMetacells <- function(seurat_obj, ...){
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- harmony::RunHarmony(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
  seurat_obj
}

#' RunUMAPMetacells
#'
#' Wrapper function to run Seurat's RunUMAP function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' NormalizeMetadata
RunUMAPMetacells <- function(seurat_obj, ...){
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::RunUMAP(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
  seurat_obj
}


#' DimPlotMetacells
#'
#' Wrapper function to run Seurat's DimPlot function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' NormalizeMetadata
DimPlotMetacells <- function(seurat_obj, ...){
  Seurat::DimPlot(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
}
