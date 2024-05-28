#' ConstructNetwork
#'
#' @description Constructs a co-expression network and groups genes into modules
#' given a Seurat object that has been prepared for network analysis.
#'
#' @return seurat_obj with the co-expression network and gene modules computed for the selected wgcna experiment
#'
#' @param seurat_obj A Seurat object
#' @param soft_power the soft power used for network construction. Automatically selected by default.
#' @param min_power the smallest soft power to be selected if soft_power=NULL
#' @param tom_outdir path to the directory where the TOM will be written
#' @param tom_name prefix name given to the TOM output file
#' @param consensus flag indicating whether or not to perform Consensus network analysis
#' @param wgcna_name name of the WGCNA experiment
#' @param ... additional parameters passed to blockwiseConsensusModules
#'
#' @details
#' ConstructNetwork builds a co-expression network and identifies clusters of
#' highly co-expressed genes (modules) from the metacell or metaspot
#' expression matrix stored in the Seurat object. Before running this function,
#' the following functions must be run on the input Seurat object:
#'
#' 1. SetupForWGCNA
#' 2. MetacellsByGroups or MetaspotsByGroups, and NormalizeMetacells
#' 3. SetDatExpr or SetMultiExpr
#' 4. TestSoftPowers or TestSoftPowersConsensus
#'
#' This function can also be used to perform consensus network analysis if consensus=TRUE
#' is selected. ConstructNetwork calls the WGCNA function blockwiseConsensusModules
#' to compute the adjacency matrix, topological overlap matrix, and to run the
#' Dynamic Tree Cut algorithm to identify gene modules. blockwiseConsensusModules
#' has numerous parameters but here we have selected default parameters that we
#' have found to provide reasonable results on a variety of single-cell and
#' spatial transcriptomic datasets.
#'
#' @import WGCNA
#' @import Seurat
#' @export
ConstructNetwork <- function(
  seurat_obj, soft_power=NULL, min_power=3,
  tom_outdir="TOM",
  tom_name = NULL,
  consensus = FALSE,
  overwrite_tom = FALSE,
  wgcna_name = NULL,
  blocks=NULL, maxBlockSize=30000, randomSeed=12345, corType="pearson",
  consensusQuantile=0.3, networkType = "signed", TOMType = "signed",
  TOMDenom = "min", scaleTOMs = TRUE, scaleQuantile = 0.8,
  sampleForScaling = TRUE, sampleForScalingFactor = 1000,
  useDiskCache = TRUE, chunkSize = NULL,
  deepSplit = 4, pamStage=FALSE, detectCutHeight = 0.995, minModuleSize = 50,
  mergeCutHeight = 0.2, saveConsensusTOMs = TRUE, ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # suffix for the tom
  if(is.null(tom_name)){
    tom_name <- gsub(' ', '_', wgcna_name)
  }

  # tom file name
  renamed <- paste0(tom_outdir, '/', tom_name, '_TOM.rda')
  if(file.exists(renamed)){
    if(overwrite_tom){
      warning(paste0('Overwriting TOM ', renamed))
    } else{
      stop(paste0("TOM ", renamed, " already exists. Set overwrite_tom = TRUE or change tom_name to proceed."))
    }
  }

  # constructing network on multiple datasets (consensus WGCNA)
  if(consensus){

    multiExpr <- GetMultiExpr(seurat_obj, wgcna_name)
    checkSets(multiExpr) # check data size

  # constructing network on a single dataset
  } else{

    # get datExpr from seurat object
    datExpr <- GetDatExpr(seurat_obj, wgcna_name)

    if(is.null(wgcna_name)){
      wgcna_name <- 'all'
    }

    nSets = 1
    setLabels = gsub(' ', '_', wgcna_name)
    shortLabels = setLabels
    multiExpr <- list()
    multiExpr[[wgcna_name]] <- list(data=datExpr)
    checkSets(multiExpr) # check data size
  }

  # make output dir for the TOM
  if(!dir.exists(tom_outdir)){
    dir.create(tom_outdir)
  }

  if(is.null(soft_power) & !consensus){
    soft_power <- GetPowerTable(seurat_obj) %>% subset(SFT.R.sq >= 0.8 & Power > min_power) %>% .$Power %>% min
    cat(paste0("Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = ", soft_power, "\n"))
  } else if(is.null(soft_power)){
    power_tables <- GetPowerTable(seurat_obj) %>% dplyr::group_split(group)
    soft_power <- sapply(power_tables, function(power_table){
      power_table %>% subset(SFT.R.sq >= 0.8 & Power > min_power) %>% .$Power %>% min
    })
    cat(paste0("Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = c(", paste0(soft_power, collapse=','), ")\n"))
  }

  # construct the network
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
    consensusTOMFilePattern = paste0(tom_outdir, '/', tom_name, "_block.%b.rda"), 
    ...
    )

  # rename consensusTOM file:
  file.rename(paste0(tom_outdir, '/', tom_name, '_block.1.rda'), renamed)

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

  # append working directory to the TOM file so it has the full path:
  net$TOMFiles <- paste0(getwd(), '/', renamed)

  # add network to seurat obj
  seurat_obj <- SetNetworkData(seurat_obj, net, wgcna_name)

  # set the modules df in the Seurat object
  mods <- GetNetworkData(seurat_obj)$colors
  seurat_obj <- SetModules(
    seurat_obj, mod_df = data.frame(
      "gene_name" = names(mods),
      "module" = factor(mods, levels=unique(mods)),
      "color" = mods
    )
  )

  seurat_obj

}
