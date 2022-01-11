
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
SelectNetworkGenes <- function(seurat_obj, gene_select="variable", fraction=0.05, gene_list=NULL, wgcna_name=NULL){

  # validate inputs:
  if(!(gene_select %in% c("variable", "fraction", "all", "custom"))){
    stop(paste0("Invalid selection gene_select: ", gene_select, '. Valid gene_selects are variable, fraction, all, or custom.'))
  }

  # handle different selection strategies
  if(gene_select == "fraction"){

    # binarize counts matrix
    expr_mat <- GetAssayData(seurat_obj, slot='counts')
    expr_mat[expr_mat>0] <- 1

    # identify genes that are expressed in at least some fraction of cells
    gene_filter <- rowSums(expr_mat) >= round(fraction*ncol(seurat_obj));
    gene_list <- rownames(seurat_obj)[gene_filter]

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
#' @keywords scRNA-seq
#' @export
#' @examples
#' SetupForWGCNA(pbmc)
SetupForWGCNA <- function(seurat_obj, wgcna_name, ...){

  # set the active WGCNA variable
  seurat_obj <- SetActiveWGCNA(seurat_obj, wgcna_name)

  # select genes for WGCNA:
  seurat_obj <- SelectNetworkGenes(seurat_obj, ...)

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
TestSoftPowers <- function(
  seurat_obj,
  use_metacells = TRUE,
  powers=c(seq(1,10,by=1), seq(12,30, by=2)),
  make_plot=TRUE, outfile="softpower.pdf", figsize=c(7,7)
){

  # add datExpr if not already added:
  if(!("datExpr" %in% names(GetActiveWGCNA(seurat_obj)))){
    seurat_obj <- SetDatExpr(seurat_obj, use_metacells)
  }

  # get datExpr
  datExpr <- GetDatExpr(seurat_obj)

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
  if(make_plot){
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
  }

  # set the power table in Seurat object:
  seurat_obj <- SetPowerTable(seurat_obj, powerTable$data)

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
  seurat_obj, soft_power=NULL, use_metacells=TRUE, cur_celltype='net', tom_outdir="TOM",
  blocks=NULL, maxBlockSize=30000, randomSeed=12345, corType="pearson",
  consensusQuantile=0.3, networkType = "signed", TOMType = "unsigned",
  TOMDenom = "min", scaleTOMs = TRUE, scaleQuantile = 0.8,
  sampleForScaling = TRUE, sampleForScalingFactor = 1000,
  useDiskCache = TRUE, chunkSize = NULL,
  deepSplit = 4, pamStage=FALSE, detectCutHeight = 0.995, minModuleSize = 50,
  mergeCutHeight = 0.2, saveConsensusTOMs = TRUE, ...
){

  # add datExpr if not already added:
  if(!("datExpr" %in% names(GetActiveWGCNA(seurat_obj)))){
    seurat_obj <- SetDatExpr(seurat_obj, use_metacells)
  }

  # get datExpr
  datExpr <- GetDatExpr(seurat_obj)

  # Add functionality to accept multiExpr and perform consensus WGCNA
  # TODO

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

  params <- list(
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

  # add parameters:
  seurat_obj <- SetWGCNAParams(seurat_obj, params)

  # add network to seurat obj
  seurat_obj <- SetNetworkData(seurat_obj, net)

  # set the modules df in the Seurat object
  mods <- GetNetworkData(seurat_obj)$colors
  seurat_obj <- SetModules(
    seurat_obj, mod_df = data.frame(
      "gene_name" = names(mods),
      "module" = mods,
      "color" = mods
    )
  )

  seurat_obj

}

#' ComputeModuleEigengene
#'
#' Internal helper function that computes module eigengene for a single module.
#'
#' @param seurat_obj A Seurat object
#' @param cur_mod name of a module found in seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_net$colors
#' @keywords scRNA-seq
#' @export
#' @examples
#' ConstructNetwork(pbmc)
ComputeModuleEigengene <- function(seurat_obj, cur_mod, group.by.vars=NULL, verbose=TRUE, wgcna_name=NULL, ...){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get module df:
  modules <- GetModules(seurat_obj, wgcna_name)

  # get genes in this module:
  cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

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
ModuleEigengenes <- function(seurat_obj, group.by.vars=NULL, verbose=TRUE, wgcna_name=NULL, ...){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # are we going to run Harmony?
  harmonized = !is.null(group.by.vars)

  me_list <- list()
  harmonized_me_list <- list()
  modules <- GetModules(seurat_obj, wgcna_name)

  # loop over modules:
  for(cur_mod in unique(modules$module)){

    print(cur_mod)

    # compute module eigengenes for this module
    seurat_obj <- ComputeModuleEigengene(
      seurat_obj = seurat_obj, cur_mod = cur_mod,
      group.by.vars=group.by.vars, verbose=verbose,
      wgcna_name, ...
    )

    # add module eigengene to ongoing list
    cur_me <- seurat_obj@reductions$ME@cell.embeddings[,1]
    me_list[[cur_mod]] <- cur_me

    # run harmony
    if(harmonized){
      # add module eigengene to ongoing list
      cur_harmonized_me <- seurat_obj@reductions$ME_harmony@cell.embeddings[,1]
      harmonized_me_list[[cur_mod]] <- cur_harmonized_me
    }
  }

  # merge module eigengene lists into a dataframe, order modules, add to Seurat obj
  me_df <- do.call(cbind, me_list)
  me_df <- WGCNA::orderMEs(me_df)
  seurat_obj <- SetMEs(seurat_obj, me_df, harmonized=FALSE, wgcna_name)

  # merge harmonized module eigengene lists into a dataframe, add to Seurat obj
  if(!is.null(group.by.vars)){
    hme_df <- do.call(cbind, harmonized_me_list)
    hme_df <- WGCNA::orderMEs(hme_df)
    seurat_obj <- SetMEs(seurat_obj, hme_df, harmonized=TRUE, wgcna_name)
  }

  # set module factor levels based on order
  MEs <- GetMEs(cur_seurat, harmonized, wgcna_name)
  modules$module <- factor(
    as.character(modules$module),
    levels=colnames(MEs)
  )
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

  # remove temp dim reductions by setting to NULL
  seurat_obj@reductions$ME <- NULL
  seurat_obj@reductions$ME_harmony <- NULL

  # return seurat object
  seurat_obj
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
ModuleConnectivity <- function(seurat_obj, harmonized=TRUE, wgcna_name=NULL, ...){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get module df, wgcna genes, and wgcna params:
  modules <- GetModules(seurat_obj, wgcna_name)
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  # datExpr for full expression dataset
  datExpr <- t(GetAssayData(
    seurat_obj,
    assay=params$metacell_assay,
    slot=params$metacell_slot
  ))[,genes_use]


  # get MEs:
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)

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
  kMEs <- cbind(modules, kMEs)
  colnames(kMEs) <- c(colnames(modules), paste0("kME_", colnames(MEs)))

  # update the modules table in the Seurat object:
  seurat_obj <- SetModules(seurat_obj, kMEs, wgcna_name)

  seurat_obj

}
