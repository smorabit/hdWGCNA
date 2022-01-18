
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

  # TODO:
  # get genes that are expressed in a fraction of cells in select groups (ie celltypes)
  # this would reduce a bias against genes that are only expressed in underrepresented
  # cell-types

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
SetupForWGCNA <- function(
  seurat_obj, wgcna_name,
  group=NULL, features = NULL,...){

  # set the active WGCNA variable
  seurat_obj <- SetActiveWGCNA(seurat_obj, wgcna_name)

  # set the active WGCNA group, this is used for actually
  # making the co-expression network
  # TODO: I never actually use WGCNAGroup anywhere??
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
  group.by=NULL, group_name=NULL,
  powers=c(seq(1,10,by=1), seq(12,30, by=2)),
  make_plot=TRUE, outfile="softpower.pdf", figsize=c(7,7)
){

  # add datExpr if not already added:
  if(!("datExpr" %in% names(GetActiveWGCNA(seurat_obj)))){
    seurat_obj <- SetDatExpr(seurat_obj, use_metacells, group.by=group.by, group_name=group_name)
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
  seurat_obj, soft_power=NULL, use_metacells=TRUE, group.by=NULL, group_name=NULL,
  tom_outdir="TOM",
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
    seurat_obj <- SetDatExpr(seurat_obj, use_metacells, group.by=group.by, group_name=group_name)
  }

  # get datExpr
  datExpr <- GetDatExpr(seurat_obj)

  if(is.null(group_name)){
    group_name <- 'all'
  }

  # Add functionality to accept multiExpr and perform consensus WGCNA
  # TODO

  nSets = 1
  setLabels = gsub(' ', '_', group_name)
  shortLabels = setLabels
  multiExpr <- list()
  multiExpr[[group_name]] <- list(data=datExpr)
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
  file.rename('ConsensusTOM-block.1.rda', paste0('TOM/', gsub(' ', '_',group_name), '_ConsensusTOM-block.1.rda'))

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
  net$TOMFiles <- paste0(getwd(), '/TOM/', group_name, '_', net$TOMFiles)

  # add network to seurat obj
  seurat_obj <- SetNetworkData(seurat_obj, net)

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
ComputeModuleEigengene <- function(
  seurat_obj, cur_mod, modules,
  group.by.vars=NULL, verbose=TRUE,
  wgcna_name=NULL, ...
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

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
ModuleEigengenes <- function(
  seurat_obj, group.by.vars=NULL,
  modules=NULL, verbose=TRUE,
  wgcna_name=NULL, ...
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # are we going to run Harmony?
  harmonized = !is.null(group.by.vars)

  me_list <- list()
  harmonized_me_list <- list()

  # get modules from Seurat object, else use provided modules
  if(is.null(modules)){
    modules <- GetModules(seurat_obj, wgcna_name)
    projected <- FALSE
  } else{
    projected <- TRUE
  }

  # get list of modules:
  mods <- levels(modules$module)

  # loop over modules:
  for(cur_mod in mods){

    print(cur_mod)

    # compute module eigengenes for this module
    seurat_obj <- ComputeModuleEigengene(
      seurat_obj = seurat_obj, cur_mod = cur_mod,
      modules=modules, group.by.vars=group.by.vars, verbose=verbose,
      wgcna_name=wgcna_name, ...
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
  if(!projected){me_df <- WGCNA::orderMEs(me_df)}
  seurat_obj <- SetMEs(seurat_obj, me_df, harmonized=FALSE, wgcna_name)

  # merge harmonized module eigengene lists into a dataframe, add to Seurat obj
  if(!is.null(group.by.vars)){
    hme_df <- do.call(cbind, harmonized_me_list)
    if(!projected){hme_df <- WGCNA::orderMEs(hme_df)}
    seurat_obj <- SetMEs(seurat_obj, hme_df, harmonized=TRUE, wgcna_name)
  }

  # set module factor levels based on order
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  print(colnames(MEs))
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

#' ModuleExprScores
#'
#' Computes module eigengenes for all WGCNA co-expression modules
#'
#' @param seurat_obj A Seurat object
#' @param method Seurat or UCell?
#' @keywords scRNA-seq
#' @export
#' @examples
#'  ModuleExprScores (pbmc)
ModuleExprScore <- function(
  seurat_obj, n_genes = 25,
   wgcna_name=NULL, method='Seurat', ...
 ){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # use all genes in the module?
  if(n_genes == "all"){
    gene_list <- lapply(mods, function(cur_mod){
      subset(modules, module == cur_mod) %>% .$gene_name
    })
  } else{
    gene_list <- lapply(mods, function(cur_mod){
      cur <- subset(modules, module == cur_mod)
      cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
        top_n(n_genes) %>% .$gene_name
    })
  }
  names(gene_list) <- mods

  # compute Module Scores with Seurat or UCell:
  if(method == "Seurat"){
    mod_scores <- Seurat::AddModuleScore(
      seurat_obj, features=gene_list, ...
    )@meta.data
  } else if(method == "UCell"){
    mod_scores <- UCell::AddModuleScore_UCell(
      seurat_obj, features=gene_list, ...
    )@meta.data
  } else{
    stop("Invalid method selection. Valid choices are Seurat, UCell")
  }

  mod_scores <- mod_scores[,(ncol(mod_scores)-length(mods)+1):ncol(mod_scores)]

  # rename module scores:
  colnames(mod_scores) <- mods

  # order cols
  col_order <- levels(modules$module)
  col_order <- col_order[col_order != 'grey']
  mod_scores <- mod_scores[,col_order]

  # add module scores to seurat object
  seurat_obj <- SetModuleScores(seurat_obj, mod_scores, wgcna_name)

  seurat_obj

}


#' AverageModuleExpr
#'
#' Computes module eigengenes for all WGCNA co-expression modules
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples
#'  AverageModuleExpr (pbmc)
AvgModuleExpr <- function(seurat_obj, n_genes = 25, wgcna_name=NULL, ...){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get wgcna params
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  # datExpr for full expression dataset
  datExpr <- GetAssayData(
    seurat_obj,
    assay=params$metacell_assay,
    slot=params$metacell_slot
  )

  # use all genes in the module?
  if(n_genes == "all"){
    gene_list <- lapply(mods, function(cur_mod){
      subset(modules, module == cur_mod) %>% .$gene_name
    })
  } else{
    gene_list <- lapply(mods, function(cur_mod){
      cur <- subset(modules, module == cur_mod)
      cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
        top_n(n_genes) %>% .$gene_name
    })
  }
  names(gene_list) <- mods

  # for each module, compute average expression of genes:
  avg_mods <- t(do.call(rbind,lapply(gene_list, function(cur_genes){colMeans(datExpr[cur_genes,])})))

  # order cols
  col_order <- levels(modules$module)
  col_order <- col_order[col_order != 'grey']
  avg_mods <- avg_mods[,col_order]

  # add avg module expression to Seurat object
  seurat_obj <- SetAvgModuleExpr(seurat_obj, avg_mods, wgcna_name)
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

  kMEs <- WGCNA::signedKME(
    datExpr,
    MEs,
    outputColumnName = "kME",
    corFnc = "bicor",
    ...
  )

  # add module color to the kMEs table
  kMEs <- cbind(modules, kMEs)
  colnames(kMEs) <- c(colnames(modules), paste0("kME_", colnames(MEs)))

  # update the modules table in the Seurat object:
  seurat_obj <- SetModules(seurat_obj, kMEs, wgcna_name)

  seurat_obj

}

#' RunEnrichr
#'
#' Computes intramodular connectivity (kME) based on module eigengenes.
#'
#' @param seurat_obj A Seurat object
#' @param dbs List of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' RunEnrichr
RunEnrichr <- function(
  seurat_obj,
  dbs = c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021'),
  max_genes = 100,
  wgcna_name=NULL, ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)

  # exclude grey
  mods <- mods[mods != 'grey']

  # run EnrichR for loop:
  combined_output <- data.frame()
  for(i in 1:length(mods)){
  	cur_mod <- mods[i]
    cur_info <- subset(modules, module == cur_mod)
    cur_info <- cur_info[,c('gene_name', paste0('kME_', cur_mod))]
    cur_genes <- top_n(cur_info, max_genes) %>% .$gene_name %>% as.character
    enriched <- enrichR::enrichr(cur_genes, dbs)

    # collapse into one db
    for(db in names(enriched)){
      cur_df <- enriched[[db]]
      if (nrow(cur_df) > 1){
        cur_df$db <- db
        cur_df$module <- cur_mod
        combined_output <- rbind(combined_output, cur_df)
      }
    }
  }

  # add GO term table to seurat object
  seurat_obj <- SetEnrichrTable(seurat_obj, combined_output, wgcna_name)
  seurat_obj

}


#' OverlapModulesDEGs
#'
#' Computes intramodular connectivity (kME) based on module eigengenes.
#'
#' @param seurat_obj A Seurat object
#' @param dbs List of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' OverlapModulesDEGs
OverlapModulesDEGs <- function(
  seurat_obj, deg_df, wgcna_name = NULL, fc_cutoff = 0.5, ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  cell_groups <- deg_df$group %>% unique

  # get modules,
  modules <- GetModules(cur_seurat)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # subset deg_df by fold change cutoff:
  if(fc_cutoff >=0 ){
    deg_df <- subset(deg_df, avg_log2FC >= fc_cutoff)
  } else{
    deg_df <- subset(deg_df, avg_log2FC <= fc_cutoff)

    # reverse the sign of the remaining fold changes
    deg_df$avg_log2FC <- -1 * deg_df$avg_log2FC
    fc_cutoff <- -1 * fc_cutoff
  }

  # size of genome based on # genes in Seurat object:
  genome.size <- nrow(seurat_obj)

  # compute overlap:
  cur_overlap <- GeneOverlap::testGeneOverlap(
    GeneOverlap::newGeneOverlap(
      cur_module_genes,
      cur_DEGs,
      genome.size=genome.size
  ))
  or <- cur_overlap@odds.ratio
  pval <- cur_overlap@pval
  jaccard <- cur_overlap@Jaccard

  # run overlaps between module gene lists and DEG lists:
  overlap_df <- do.call(rbind, lapply(mods, function(cur_mod){
    cur_module_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name
    cur_overlap_df <- do.call(rbind, lapply(cell_groups, function(cur_group){
      cur_DEGs <- deg_df %>% subset(group == cur_group & p_val_adj <= 0.05 & avg_log2FC > fc_cutoff) %>% .$gene
      cur_overlap <- testGeneOverlap(newGeneOverlap(
          cur_module_genes,
          cur_DEGs,
          genome.size=genome.size
      ))
      c(cur_overlap@odds.ratio, cur_overlap@pval, cur_overlap@Jaccard)
    })) %>% as.data.frame
    colnames(cur_overlap_df) <- c('odds_ratio', 'pval', 'Jaccard')
    cur_overlap_df$module <- cur_mod
    cur_overlap_df$group <- cell_groups

    # module color:
    cur_overlap_df$color <- modules %>% subset(module == cur_mod) %>% .$color %>% unique
    cur_overlap_df
  }))

  # adjust for multiple comparisons:
  overlap_df$fdr <- p.adjust(overlap_df$pval, method='fdr')

  # significance level:
  overlap_df$Significance <- gtools::stars.pval(overlap_df$fdr)
  overlap_df$Significance <- ifelse(
    overlap_df$Significance == '.', '',
    overlap_df$Significance
  )

  # set factor levels for modules:
  overlap_df$module <- factor(overlap_df$module, levels=mods)

  # re-arrange columns:
  overlap_df <- overlap_df %>% select(c(module, group, color, odds_ratio, pval, fdr, Significance, Jaccard))

  overlap_df
}


#' ProjectModules
#'
#' Computes intramodular connectivity (kME) based on module eigengenes.
#'
#' @param seurat_obj A Seurat object
#' @param dbs List of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name
#' @keywords scRNA-seq
#' @export
#' @examples
#' ProjectModules
ProjectModules <- function(
  seurat_obj, seurat_ref,
  group.by.vars=NULL,
  gene_mapping=NULL, # table mapping genes from species 1 to species 2
  genome1_col=NULL, genome2_col=NULL,
  scale_genes = FALSE,
  wgcna_name=NULL, wgcna_name_proj=NULL,
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_ref@misc$active_wgcna}

  # get modules to be projected:
  modules <- GetModules(seurat_ref, wgcna_name)

  # are we mapping gene names?
  if(!is.null(gene_mapping)){
    modules <- TransferModuleGenome(modules, gene_mapping, genome1_col, genome2_col)
  }

  print('here')

  # get genes that overlap between WGCNA genes & seurat_obj genes:
  gene_names <- modules$gene_name
  genes_use <- intersect(gene_names, rownames(seurat_obj))

  # subset modules by genes in this seurat object:
  modules <- modules %>% subset(gene_name %in% genes_use)

  # setup new seurat obj for wgcna:
  if(!(wgcna_name_proj %in% names(seurat_obj@misc))){
    seurat_obj <- SetupForWGCNA(
      seurat_obj,
      wgcna_name = wgcna_name_proj,
      features = genes_use
    )
  }

  # scale the dataset if needed:
  if(!scale_genes & sum(genes_use %in% rownames(GetAssayData(seurat_obj, slot='scale.data'))) == length(genes_use)){
    print("Scaling already done.")
  } else if(scale_genes){
    print("Scaling dataset...")
    seurat_obj <- Seurat::ScaleData(
      seurat_obj, features = genes_use, ...
    )
  }


  # project modules:
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars=group.by.vars,
    modules = modules,
    wgcna_name = wgcna_name_proj,
    ...
  )

  seurat_obj

}



#' TransferModuleGenome
#'
#' Takes a module table and a gene mapping table (like from biomart) with gene
#' names from two genomes in order to switch
#'
#' @param modules module table from GetModules function
#' @param gene_mapping table that has the gene names for both genomes
#' @param genome1_col the column in gene_mapping that has gene names for the genome currently present in modules
#' @param genome2_col the column in gene_mapping that has gene names for the genome to transfer to
#' @keywords scRNA-seq
#' @export
#' @examples
#' TransferModuleGenome
TransferModuleGenome <- function(
  modules, gene_mapping,
  genome1_col=NULL, genome2_col=NULL
){

  # use the first & second columns if these are null
  if(is.null(genome1_col)){genome1_col <- colnames(gene_mapping)[1]}
  if(is.null(genome2_col)){genome2_col <- colnames(gene_mapping)[2]}

  # switch gene names to human:
  gene_mapping <- gene_mapping[,c(genome1_col, genome2_col)]

  # only keep genome1 genes that are in the WGCNA gene list:
  gene_list <- modules$gene_name
  gene_mapping <- gene_mapping[gene_mapping[,genome1_col] %in% gene_list,]

  # subset modules
  modules <- subset(modules, gene_name %in% gene_mapping[,genome1_col])

  print('here')

  # match the order of the given gene list, and remove NA entries
  gene_match <- match(gene_list, gene_mapping[,genome1_col])
  gene_mapping <- na.omit(gene_mapping[gene_match,])

  # update modules table with the new gene names
  # modules <- na.omit(modules[gene_match,])

  print('here')
  print(dim(modules))
  print(dim(gene_mapping))
  print(length(gene_match))

  modules$gene_name <- gene_mapping[,genome2_col]

  modules

}
