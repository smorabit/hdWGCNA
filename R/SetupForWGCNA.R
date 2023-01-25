
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
#' SelectNetworkGenes(pbmc)
SelectNetworkGenes <- function(
  seurat_obj,
  gene_select="variable", fraction=0.05,
  group.by=NULL, # should be a column in the Seurat object, eg clusters
  gene_list=NULL,
  wgcna_name=NULL
 ){

  # validate inputs:
  if(!(gene_select %in% c("variable", "fraction", "all", "custom"))){
    stop(paste0("Invalid selection gene_select: ", gene_select, '. Valid gene_selects are variable, fraction, all, or custom.'))
  }

  # get assay
  assay <- DefaultAssay(seurat_obj)
  tryCatch(
   tmp <- dim(seurat_obj[[assay]]@counts),
   error = function(e){
     print("Assay must contain counts slot.")
   })

  # handle different selection strategies
  if(gene_select == "fraction"){

    # binarize counts matrix in chunks to save memory
    expr_mat <- GetAssayData(seurat_obj, slot='counts')
    n_chunks <- ceiling(ncol(expr_mat) / 10000)

    if(n_chunks == 1){
      chunks <- factor(rep(1), levels=1)
    } else{
      chunks <- cut(1:nrow(expr_mat), n_chunks)
    }
    expr_mat <- do.call(rbind, lapply(levels(chunks), function(x){
      cur <- expr_mat[chunks == x,]
      cur[cur > 0] <- 1
      cur
    }))

    group_gene_list <- list()
    if(!is.null(group.by)){
      # identify genes that are expressed
      groups <- unique(seurat_obj@meta.data[[group.by]])
      for(cur_group in groups){

        # subset expression matrix by this group
        cur_expr <- expr_mat[,seurat_obj@meta.data[[group.by]] == cur_group]
        print(dim(cur_expr))

        gene_filter <- Matrix::rowSums(cur_expr) >= round(fraction*ncol(cur_expr));
        group_gene_list[[cur_group]] <- rownames(seurat_obj)[gene_filter]
      }
      gene_list <- unique(unlist(group_gene_list))

    } else{
      # identify genes that are expressed in at least some fraction of cells
      gene_filter <- Matrix::rowSums(expr_mat) >= round(fraction*ncol(seurat_obj));
      gene_list <- rownames(seurat_obj)[gene_filter]
    }

  } else if(gene_select == "variable"){
    gene_list <- VariableFeatures(seurat_obj)

  } else if(gene_select == "all"){
    gene_list <- rownames(seurat_obj)

  } else if(gene_select == "custom"){

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
  seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)

  # return updated seurat obj
  seurat_obj

}




#' SetupForWGCNA
#'
#' This function gets the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name name of the WGCNA experiment
#' @param features list of features to use for WGCNA
#' @param metacell_location name of the WGCNA experiment to copy the metacell object from
#' @param ... additional parameters to pass to SelectNetworkGenes
#' @param group (parameter not used anymore, can ignore)
#' @keywords scRNA-seq
#' @export
#' @examples
#' SetupForWGCNA(pbmc)
SetupForWGCNA <- function(
  seurat_obj, wgcna_name,
  features = NULL,
  metacell_location = NULL,
  group=NULL,
  ...
){

  # set the active WGCNA variable
  seurat_obj <- SetActiveWGCNA(seurat_obj, wgcna_name)

  # set the active WGCNA group, this is used for actually
  # making the co-expression network
  if(is.null(group)){
    seurat_obj <- SetWGCNAGroup(seurat_obj, group='all', wgcna_name)
  } else{
    seurat_obj <- SetWGCNAGroup(seurat_obj, group, wgcna_name)
  }

  # select genes for WGCNA, else use pre-selected:
  if(is.null(features)){
    seurat_obj <- SelectNetworkGenes(seurat_obj, wgcna_name=wgcna_name, ...)
  } else{
    seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list=features, wgcna_name=wgcna_name)
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
