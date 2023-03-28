
#' SelectNetworkGenes
#'
#' Select genes that will be used for co-expression network analysis
#' 
#' @param seurat_obj A Seurat object
#' @param type How to select genes? Select "variable", "fraction", "all", or "custom".
#' @param fraction A numeric that determines the minimum cells that a gene must be expressed in order to be included. For example, fraction = 0.05 means that 5% of cells must express a gene (count > 0) for it to be included.
#' @param gene_list A character string of gene names, only used if type = "custom"
#' 
#' @details 
#' SelectNetworkGenes allows us to specify the genes that will be used for co-expression network analysis.
#' This function is called by SetupForWGCNA. By default, the variable features in VariableFeatures(seurat_obj) are used.
#' A custom gene list can also be used with the gene_list parameter and setting gene_select="custom".
#' 
#' We can also identify genes that are expressed above 0 in a certain proportion of the dataset by settting gene_select='fraction'.
#' For example, by setting fraction=0.1 and group.by='seurat_clusters', this function will identify the set of genes 
#' that are expressed in 10% of cells in at least one of the clusters.
#' 
#' 
#' @export
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

    # make sure that there aren't duplicates
    gene_list <- unique(gene_list)

   # all selected features should be present in the Seurat object:
    if(!all(gene_list %in% rownames(seurat_obj))){
      stop("Some selected features are not found in rownames(seurat_obj).")
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
