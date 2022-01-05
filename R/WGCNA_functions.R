
#' SelectNetworkGenes
#'
#' This function
#' on neighboring cells in provided groupings, such as cluster or cell type.
#' @param seurat_obj A Seurat object
#' @param type How to select genes? Select "variable", "fraction", "all", or "custom".
#' @param fraction A numeric that determines the minimum cells that a gene must be expressed in order to be included. For example, fraction = 0.05 means that 5% of cells must express a gene (count > 0) for it to be included.
#' @param gene_list A character string of gene names, only used if type = "custom"
#' @keywords scRNA-seq
#' @export
#' @examples
#' MetacellsByGroups(pbmc)
SelectNetworkGenes <- function(seurat_obj, type="variable", fraction=0.05, gene_list=NULL){

  # validate inputs:
  if(!(type %in% c("variable", "fraction", "all", "custom"))){
    stop(paste0("Invalid selection type: ", type, '. Valid types are variable, fraction, all, or custom.'))
  }

  # handle different selection strategies
  if(type == "fraction"){

    # binarize counts matrix
    expr_mat <- GetAssayData(cur_seurat, slot='counts')
    expr_mat[expr_mat>0] <- 1

    # identify genes that are expressed in at least some fraction of cells
    gene_filter <- rowSums(expr_mat) >= round(fraction*ncol(cur_seurat));
    gene_list <- rownames(seurat_obj)[gene_filter]

  } else if(type == "variable"){
    gene_list <- VariableFeatures(seurat_obj)

  } else if(type == "all"){
    gene_list <- rownames(seurat_obj)

  } else if(type == "custom"){

    gene_list <- gene_list

    # check that custom genes are present in the seurat object:
    check_genes <- gene_list %in% rownames(seurat_obj)
    if(sum(check_genes) < length(gene_list)){
      stop(paste("Some genes not present in seurat object:", paste(gene_list[!check_genes], collapse=', ')))
    }
  }

  # make sure there's more than 0 genes:
  if(length(gene_list) == 0){
    stop("No genes found")
  }

  # throw a warning if there's very few genes:
  if(length(gene_list) <= 100){
    warning(paste0("Very few genes selected (", length(gene_list), "), perhaps use a different method to select genes."))
  }

  # add genes to gene list
  seurat_obj@misc$wgcna_genes <- gene_list

  # set VariableFeatures slot for metacell object as these genes
  VariableFeatures(seurat_obj@misc$wgcna_metacell_obj) <- gene_list

  # return updated seurat obj
  seurat_obj

}
