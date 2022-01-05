
#' NormalizeMetacells
#'
#' Wrapper function to run Seurat's NormalizeData function on the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' NormalizeMetadata
NormalizeMetacells <- function(seurat_obj, ...){
  seurat_obj@misc$wgcna_metacell_obj <- Seurat::NormalizeData(seurat_obj@misc$wgcna_metacell_obj, ...)
  seurat_obj
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
ScaleMetacells <- function(seurat_obj, ...){
  if(!exists("features")){
    features = VariableFeatures(seurat_obj)
  }
  seurat_obj@misc$wgcna_metacell_obj <- Seurat::ScaleData(seurat_obj@misc$wgcna_metacell_obj, ...)
  seurat_obj
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
RunPCAMetacells <- function(seurat_obj, ...){
  seurat_obj@misc$wgcna_metacell_obj <- Seurat::RunPCA(seurat_obj@misc$wgcna_metacell_obj, ...)
  seurat_obj
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
  seurat_obj@misc$wgcna_metacell_obj <- harmony::RunHarmony(seurat_obj@misc$wgcna_metacell_obj, ...)
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
  seurat_obj@misc$wgcna_metacell_obj <- Seurat::RunUMAP(seurat_obj@misc$wgcna_metacell_obj, ...)
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
  Seurat::DimPlot(seurat_obj@misc$wgcna_metacell_obj, ...)
}
