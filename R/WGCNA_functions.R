
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

#' GetDatExpr
#'
#' This function gets the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' GetDatExpr(pbmc)
GetDatExpr <- function(seurat_obj, ...){
  params <- seurat_obj@misc$wgcna_params
  genes_use <- seurat_obj@misc$wgcna_genes
  assay <- params$metacell_assay
  datExpr <- as.data.frame(
    Seurat::GetAssayData(
      seurat_obj,
      assay=assay,
      slot='data'
    )[genes_use,]
  )
  tic("transposing dataset")
  datExpr <- as.data.frame(t(datExpr))
  toc()
  # only get good genes:
  tic("Running goodGenes")
  datExpr <- datExpr[,WGCNA::goodGenes(datExpr, ...)]
  toc()

  datExpr
}


#' SetupForWGCNA
#'
#' This function gets the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' SetupForWGCNA(pbmc)
SetupForWGCNA <- function(seurat_obj){

  # setup datExpr for WGCNA
  seurat_obj@misc$wgcna_datExpr <- GetDatExpr(seurat_obj)
  seurat_obj
}

#' TestSoftPowers
#'
#' This function gets the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @param powers
#' @param outfile filepath for output pdf to be generated. Default is
#' @param figsize numeric determining the height and width of the output figure. Default is c(7,7)
#' @keywords scRNA-seq
#' @export
#' @examples
#' TestSoftPowers(pbmc)
TestSoftPowers <- function(seurat_obj, powers=c(seq(1,10,by=1), seq(12,30, by=2)), outfile="softpower.pdf", figsize=c(7,7)){

  # get datExpr from seurat object
  datExpr <- seurat_obj@misc$wgcna_datExpr

  # Call the network topology analysis function for each set in turn
  powerTable = list(
    data = WGCNA::pickSoftThreshold(
      datExpr,
      powerVector=powers,
      verbose = 100,
      networkType="signed",
      corFnc="bicor"
    )[[2]]
  );

  # Plot the results:
  pdf(outfile, height=figsize[1], width=figsize[2], useDingbats=FALSE)

      colors = c("blue", "red","black")
      # Will plot these columns of the returned scale free analysis tables
      plotCols = c(2,5,6,7)
      colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Mean connectivity",
      "Max connectivity");

      # Get the minima and maxima of the plotted points
      ylim = matrix(NA, nrow = 2, ncol = 4);
      for (col in 1:length(plotCols)){
        ylim[1, col] = min(ylim[1, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
        ylim[2, col] = max(ylim[2, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
      }

      # Plot the quantities in the chosen columns vs. the soft thresholding power
      par(mfcol = c(2,2));
      par(mar = c(4.2, 4.2 , 2.2, 0.5))
      cex1 = 0.7;

      for (col in 1:length(plotCols)){
        plot(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
        xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
        main = colNames[col]);
        addGrid();

        if (col==1){
          text(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
          labels=powers,cex=cex1,col=colors[1]);
        } else
        text(powerTable$data[,1], powerTable$data[,plotCols[col]],
        labels=powers,cex=cex1,col=colors[1]);
      }
  dev.off()

  seurat_obj@misc$wgcna_powerTable <- powerTable$data
  seurat_obj
}

#' ConstructNetwork
#'
#' This function gets the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @param power
#' @keywords scRNA-seq
#' @export
#' @examples
#' ConstructNetwork(pbmc)
ConstructNetwork <- function(seurat_obj, soft_power=NULL, cur_celltype='net'){

  # get datExpr from seurat object
  datExpr <- seurat_obj@misc$wgcna_datExpr

  nSets = 1
  setLabels = gsub(' ', '_', cur_celltype)
  shortLabels = setLabels
  multiExpr <- list()
  multiExpr[[cur_celltype]] <- list(data=datExpr)
  checkSets(multiExpr) # check data size



  net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                                           maxBlockSize = 30000, ## This should be set to a smaller size if the user has limited RAM
                                           randomSeed = 12345,
                                           corType = "pearson", ## no use for bicor
                                           power = soft_power,
                                           consensusQuantile = 0.3,
                                           networkType = "signed",
                                           TOMType = "unsigned",
                                           TOMDenom = "min",
                                           scaleTOMs = TRUE, scaleQuantile = 0.8,
                                           sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                                           useDiskCache = TRUE, chunkSize = NULL,
                                           deepSplit = 4,
                                           pamStage=FALSE,
                                           detectCutHeight = 0.995, minModuleSize = 50,
                                           mergeCutHeight = 0.2,
                                           saveConsensusTOMs = TRUE,
                                           consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")


  # rename consensusTOM file:
  file.rename('ConsensusTOM-block.1.rda', paste0('data/', gsub(' ', '_',cur_celltype), '_ConsensusTOM-block.1.rda'))

  seurat_obj@misc$wgcna_net <- net
  seurat_obj

}
