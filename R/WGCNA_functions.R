
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
    expr_mat <- GetAssayData(seurat_obj, slot='counts')
    expr_mat[expr_mat>0] <- 1

    # identify genes that are expressed in at least some fraction of cells
    gene_filter <- rowSums(expr_mat) >= round(fraction*ncol(seurat_obj));
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

  # get parameters from seurat object
  params <- seurat_obj@misc$wgcna_params
  genes_use <- seurat_obj@misc$wgcna_genes
  assay <- params$metacell_assay

  # get expression data from seurat obj
  datExpr <- as.data.frame(
    Seurat::GetAssayData(
      seurat_obj,
      assay=assay,
      slot='data'
    )[genes_use,]
  )

  # transpose data
  datExpr <- as.data.frame(t(datExpr))

  # only get good genes:
  datExpr <- datExpr[,WGCNA::goodGenes(datExpr, ...)]
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
#' This function constructs a co-expression network from a Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param soft_power
#' @keywords scRNA-seq
#' @export
#' @examples
#' ConstructNetwork(pbmc)
ConstructNetwork <- function(
  seurat_obj, soft_power=NULL, cur_celltype='net', tom_outdir="TOM",
  blocks=NULL, maxBlockSize=30000, randomSeed=12345, corType="pearson",
  consensusQuantile=0.3, networkType = "signed", TOMType = "unsigned",
  TOMDenom = "min", scaleTOMs = TRUE, scaleQuantile = 0.8,
  sampleForScaling = TRUE, sampleForScalingFactor = 1000,
  useDiskCache = TRUE, chunkSize = NULL,
  deepSplit = 4, pamStage=FALSE, detectCutHeight = 0.995, minModuleSize = 50,
  mergeCutHeight = 0.2, saveConsensusTOMs = TRUE, ...
){

  # get datExpr from seurat object
  datExpr <- seurat_obj@misc$wgcna_datExpr

  nSets = 1
  setLabels = gsub(' ', '_', cur_celltype)
  shortLabels = setLabels
  multiExpr <- list()
  multiExpr[[cur_celltype]] <- list(data=datExpr)
  checkSets(multiExpr) # check data size

  if(!dir.exists(tom_outdir)){
    dir.create(tom_outdir)
  }


  net <- WGCNA::blockwiseConsensusModules(
    multiExpr,
    power = soft_power,
    blocks = blocks,
    maxBlockSize = maxBlockSize, ## This should be set to a smaller size if the user has limited RAM
    randomSeed = randomSeed,
    corType = corType,
    consensusQuantile = consensusQuantile,
    networkType = networkType,
    TOMType = TOMType,
    TOMDenom = TOMDenom,
    scaleTOMs = scaleTOMs, scaleQuantile = scaleQuantile,
    sampleForScaling = sampleForScaling, sampleForScalingFactor = sampleForScalingFactor,
    useDiskCache = useDiskCache, chunkSize = chunkSize,
    deepSplit = deepSplit,
    pamStage=pamStage,
    detectCutHeight = detectCutHeight, minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    saveConsensusTOMs = saveConsensusTOMs,
    consensusTOMFilePattern = "ConsensusTOM-block.%b.rda", ...)

  # rename consensusTOM file:
  file.rename('ConsensusTOM-block.1.rda', paste0('TOM/', gsub(' ', '_',cur_celltype), '_ConsensusTOM-block.1.rda'))

  # add network parameters to the Seurat object:

  net_params <- list(
    power = soft_power,
    blocks = blocks,
    maxBlockSize = maxBlockSize, ## This should be set to a smaller size if the user has limited RAM
    randomSeed = randomSeed,
    corType = corType,
    consensusQuantile = consensusQuantile,
    networkType = networkType,
    TOMType = TOMType,
    TOMDenom = TOMDenom,
    scaleTOMs = scaleTOMs, scaleQuantile = scaleQuantile,
    sampleForScaling = sampleForScaling, sampleForScalingFactor = sampleForScalingFactor,
    useDiskCache = useDiskCache, chunkSize = chunkSize,
    deepSplit = deepSplit,
    pamStage=pamStage,
    detectCutHeight = detectCutHeight, minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    saveConsensusTOMs = saveConsensusTOMs
  )

  seurat_obj@misc$wgcna_params$net_params <- net_params

  # add network to seurat obj
  seurat_obj@misc$wgcna_net <- net
  seurat_obj

}

#' ComputeModuleEigengene
#'
#' Internal helper function that computes module eigengene for a single module.
#'
#' @param seurat_obj A Seurat object
#' @param cur_mod name of a module found in seurat_obj@misc$wgcna_net$colors
#' @keywords scRNA-seq
#' @export
#' @examples
#' ConstructNetwork(pbmc)
ComputeModuleEigengene <- function(seurat_obj, cur_mod, group.by.vars=NULL, verbose=TRUE, ...){

  # get genes in this module
  cur_genes <- names(seurat_obj@misc$wgcna_net$colors[seurat_obj@misc$wgcna_net$colors == cur_mod])

  # run PCA with Seurat function
  cur_pca <- Seurat::RunPCA(
    seurat_obj,
    features = cur_genes,
    reduction.key=paste0('pca', cur_mod),
    verbose=verbose, ...
  )@reductions$pca@cell.embeddings

  # add this PCA as its own reduction in the seurat object
  seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
    embeddings = cur_pca,
    assay = Seurat::DefaultAssay(seurat_obj)
  )

  # run harmony
  if(!is.null(group.by.vars)){
    cur_harmony <- harmony::RunHarmony(
      seurat_obj,
      group.by.vars=group.by.vars,
      reduction="ME", verbose=verbose, ...
    )@reductions$harmony@cell.embeddings

    # add harmonized PCA as its own reduction in the seurat object
    seurat_obj@reductions$ME_harmony <- Seurat::CreateDimReducObject(
      embeddings = cur_harmony,
      assay = Seurat::DefaultAssay(seurat_obj)
    )
  }

  # return seurat object
  seurat_obj

}

#' ModuleEigengenes
#'
#' Computes module eigengenes for all WGCNA co-expression modules
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#'  ModuleEigengenes(pbmc)
ModuleEigengenes <- function(seurat_obj, group.by.vars=NULL, verbose=TRUE, ...){

  me_list <- list()
  harmonized_me_list <- list()
  modules <- unique(seurat_obj@misc$wgcna_net$colors)

  # loop over modules:
  for(cur_mod in modules){

    print(cur_mod)

    # compute module eigengenes for this module
    seurat_obj <- ComputeModuleEigengene(
      seurat_obj = seurat_obj, cur_mod = cur_mod,
      group.by.vars=group.by.vars, verbose=verbose, ...
    )

    # add module eigengene to ongoing list
    cur_me <- seurat_obj@reductions$ME@cell.embeddings[,1]
    me_list[[cur_mod]] <- cur_me

    # run harmony
    if(!is.null(group.by.vars)){
      # add module eigengene to ongoing list
      cur_harmonized_me <- seurat_obj@reductions$ME_harmony@cell.embeddings[,1]
      harmonized_me_list[[cur_mod]] <- cur_harmonized_me
    }
  }

  # merge module eigengene lists into a dataframe, add to Seurat obj
  me_df <- do.call(cbind, me_list)
  seurat_obj@misc$wgcna_MEs <- me_df

  # merge harmonized module eigengene lists into a dataframe, add to Seurat obj
  if(!is.null(group.by.vars)){
    hme_df <- do.call(cbind, harmonized_me_list)
    seurat_obj@misc$wgcna_hMEs <- hme_df
  }

  # remove temp dim reductions by setting to NULL
  seurat_obj@reductions$ME <- NULL
  seurat_obj@reductions$ME_harmony <- NULL

  # return seurat object
  seurat_obj
}

#' GetMEs
#'
#' Function to retrieve module eigengens from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#'  ModuleEigengenes(pbmc)
GetMEs <- function(seurat_obj, harmonized=TRUE){
  if(harmonized == TRUE && !is.null(seurat_obj@misc$wgcna_hMEs)){
    MEs <- seurat_obj@misc$wgcna_hMEs
  } else{
    MEs <- seurat_obj@misc$wgcna_MEs
  }
  MEs
}


#' ModuleConnectivity
#'
#' Computes intramodular connectivity (kME) based on module eigengenes.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleConnectivity(pbmc)
ModuleConnectivity <- function(seurat_obj, harmonized=TRUE, ...){

  # get expression matrix
  genes_use <- names(seurat_obj@misc$wgcna_net$colors)

  # datExpr for full expression dataset
  # transpose expression data and subset by genes used for WGCNA
  datExpr <- t(GetAssayData(
    seurat_obj,
    assay=seurat_obj@misc$wgcna_params$metacell_assay,
    slot=seurat_obj@misc$wgcna_params$metacell_slot
  ))[,genes_use]

  # get MEs:
  MEs <- GetMEs(seurat_obj=seurat_obj, harmonized=harmonized)

  tic("SignedKME")
  kMEs <- signedKME(
    datExpr,
    MEs,
    outputColumnName = "kME",
    corFnc = "bicor",
    ...
  )
  toc()

  # add module color to the kMEs table
  kMEs <- cbind(cur_seurat@misc$wgcna_net$colors, kMEs)
  colnames(kMEs) <- c('module', colnames(MEs))

  seurat_obj@misc$wgcna_kMEs <- kMEs
  seurat_obj

}
