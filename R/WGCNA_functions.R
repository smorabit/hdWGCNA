
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

#' TestSoftPowers
#'
#' Compute the scale-free topology model fit for different soft power thresholds
#'
#' @param seurat_obj A Seurat object
#' @param powers numeric vector specifying soft powers to test
#' @param use_metacells logical flag for whether to use the metacell expression matrix
#' @param networkType The type of network to use for network analysis. Options are "signed" (default), "unsigned", or "signed hybrid". This should be consistent with the network chosen for ConstructNetwork
#' @param corFnc Correlation function for the gene-gene correlation adjacency matrix.
#' @param setDatExpr logical flag indicating whether to run setDatExpr.
#' @param group.by A string containing the name of a column in the Seurat object with cell groups (clusters, cell types, etc). If NULL (default), hdWGCNA uses the Seurat Idents as the group.
#' @param group_name A string containing a group present in the provided group.by column or in the Seurat Idents. A character vector can be provided to select multiple groups at a time.
#' @param ... additional parameters passed to SetDatExpr
#' @keywords scRNA-seq
#' @export
#' @examples
#' TestSoftPowers(pbmc)
TestSoftPowers <- function(
  seurat_obj,
  powers=c(seq(1,10,by=1), seq(12,30, by=2)),
  use_metacells = TRUE,
  networkType="signed",
  corFnc='bicor',
  setDatExpr = FALSE,
  group.by=NULL, group_name=NULL,
  ...
){

  # add datExpr if not already added:
  if(!("datExpr" %in% names(GetActiveWGCNA(seurat_obj))) | setDatExpr == TRUE){
    seurat_obj <- SetDatExpr(seurat_obj, use_metacells=use_metacells, group.by=group.by, group_name=group_name)
  }

  # get datExpr
  datExpr <- GetDatExpr(seurat_obj)

  # Call the network topology analysis function for each set in turn
  powerTable = list(
    data = WGCNA::pickSoftThreshold(
      datExpr,
      powerVector=powers,
      verbose = 100,
      networkType=networkType,
      corFnc=corFnc,
      ...
    )[[2]]
  );

  # set the power table in Seurat object:
  seurat_obj <- SetPowerTable(seurat_obj, powerTable$data)

  seurat_obj
}



#' TestSoftPowersConsensus
#'
#' Compute the scale-free topology model fit for different soft power thresholds separately for each input dataset
#'
#' @param seurat_obj A Seurat object
#' @param powers numeric vector specifying soft powers to test
#' @param use_metacells logical flag for whether to use the metacell expression matrix
#' @param networkType The type of network to use for network analysis. Options are "signed" (default), "unsigned", or "signed hybrid". This should be consistent with the network chosen for ConstructNetwork
#' @param corFnc Correlation function for the gene-gene correlation adjacency matrix.
#' @param setDatExpr logical flag indicating whether to run setDatExpr.
#' @param group.by A string containing the name of a column in the Seurat object with cell groups (clusters, cell types, etc). If NULL (default), hdWGCNA uses the Seurat Idents as the group.
#' @param group_name A string containing a group present in the provided group.by column or in the Seurat Idents. A character vector can be provided to select multiple groups at a time.
#' @param multi.group.by A string containing the name of a column in the Seurat object with groups for consensus WGCNA (dataset, sample, condition, etc)
#' @param multi_groups A character vecrtor containing the names of groups to select
#' @param ... additional parameters passed to SetDatExpr
#' @keywords scRNA-seq
#' @export
#' @examples
#' # TestSoftPowers(pbmc)
TestSoftPowersConsensus <- function(
  seurat_obj,
  powers=c(seq(1,10,by=1), seq(12,30, by=2)),
  use_metacells = TRUE,
  networkType="signed",
  corFnc='bicor',
  setDatExpr = FALSE,
  group.by=NULL, group_name=NULL,
  multi.group.by = NULL,
  multi_groups = NULL,
  ...
){

  # add multiExpr if not already added:
  if(!("multiExpr" %in% names(GetActiveWGCNA(seurat_obj))) | setDatExpr == TRUE){

    # set
    seurat_obj <- SetMultiExpr(
      seurat_obj,
      group_name = group_name,
      group.by = group.by,
      multi.group.by = multi.group.by,
      multi_groups = multi_groups
    )
  }

  multiExpr <- GetMultiExpr(seurat_obj)

  # pick soft thresh for each consensus group:
  powerTables <- list()
  for(i in 1:length(multiExpr)){

    cur_group <- names(multiExpr)[i]
    print(cur_group)

    # Call the network topology analysis function for each set in turn
    powerTable = list(
      data = WGCNA::pickSoftThreshold(
        multiExpr[[cur_group]]$data,
        powerVector=powers,
        verbose = 100,
        networkType=networkType,
        corFnc=corFnc,
        ...
      )[[2]]
    );
    powerTable$data$group <- cur_group
    powerTables[[cur_group]] <- powerTable$data

  }

  # merge the power tables
  powerTable <- do.call(rbind, powerTables)

  # set the power table in Seurat object:
  seurat_obj <- SetPowerTable(seurat_obj, powerTable)
  seurat_obj
}



#' ConstructNetwork
#'
#' This function constructs a co-expression network from a Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param soft_power the soft power used for network construction. Automatically selected by default.
#' @param tom_outdir path to the directory where the TOM will be written
#' @param ... additional parameters passed to SetDatExpr and blockwiseConsensusModules
#' @keywords scRNA-seq
#' @export
#' @examples
#' ConstructNetwork(pbmc)
ConstructNetwork <- function(
  seurat_obj, soft_power=NULL, min_power=3,
  tom_outdir="TOM",
  use_metacells=TRUE,
  setDatExpr=FALSE, group.by=NULL,
  group_name=NULL,
  tom_name = NULL,
  consensus = FALSE,
  multi.group.by = NULL,
  multi_groups = NULL,
  overwrite_tom = FALSE,
  blocks=NULL, maxBlockSize=30000, randomSeed=12345, corType="pearson",
  consensusQuantile=0.3, networkType = "signed", TOMType = "signed",
  TOMDenom = "min", scaleTOMs = TRUE, scaleQuantile = 0.8,
  sampleForScaling = TRUE, sampleForScalingFactor = 1000,
  useDiskCache = TRUE, chunkSize = NULL,
  deepSplit = 4, pamStage=FALSE, detectCutHeight = 0.995, minModuleSize = 50,
  mergeCutHeight = 0.2, saveConsensusTOMs = TRUE, ...
){

  wgcna_name <- seurat_obj@misc$active_wgcna


  # suffix for the tom
  if(is.null(tom_name)){
    tom_name <- gsub(' ', '_', group_name)
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

    # add multiExpr if not already added:
    if(!("multiExpr" %in% names(GetActiveWGCNA(seurat_obj))) | setDatExpr == TRUE){

      # set
      seurat_obj <- SetMultiExpr(
        seurat_obj,
        group_name = group_name,
        group.by = group.by,
        multi.group.by = multi.group.by,
        multi_groups = multi_groups
      )
    }

    multiExpr <- GetMultiExpr(seurat_obj)
    checkSets(multiExpr) # check data size

  # constructing network on a single dataset
  } else{

    # add datExpr if not already added:
    if(!("datExpr" %in% names(GetActiveWGCNA(seurat_obj))) | setDatExpr == TRUE){
      seurat_obj <- SetDatExpr(
        seurat_obj,
        group_name = group_name,
        group.by=group.by,
        use_metacells=use_metacells,
        return_seurat=TRUE
       )

    }

    # get datExpr from seurat object
    datExpr <- GetDatExpr(seurat_obj)

    if(is.null(group_name)){
      group_name <- 'all'
    }
    #
    # # suffix for the tom
    # if(is.null(tom_name)){
    #   tom_name <- gsub(' ', '_', group_name)
    # }
    #
    # # tom file name
    # renamed <- paste0(tom_outdir, '/', tom_name, '_TOM.rda')
    # if(file.exists(renamed)){
    #   if(overwrite_tom){
    #     warning(paste0('Overwriting TOM ', renamed))
    #   } else{
    #     stop(paste0("TOM ", renamed, " already exists. Set overwrite_tom = TRUE or change tom_name to proceed."))
    #   }
    # }

    nSets = 1
    setLabels = gsub(' ', '_', group_name)
    shortLabels = setLabels
    multiExpr <- list()
    multiExpr[[group_name]] <- list(data=datExpr)
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
    consensusTOMFilePattern = "ConsensusTOM-block.%b.rda", ...)

  # rename consensusTOM file:
  file.rename('ConsensusTOM-block.1.rda', renamed)

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


#' ComputeModuleEigengene
#'
#' Internal helper function that computes module eigengene for a single module.
#'
#' @param seurat_obj A Seurat object
#' @param cur_mod name of a module found in seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_net$colors
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are "linear", "poisson", or "negbinom"
#' @param pc_dim Which PC to use as the module eigengene? Default to 1.
#' @param assay Assay in seurat_obj to compute module eigengenes. Default is DefaultAssay(seurat_obj)
#' @param wgcna_name name of the WGCNA experiment
#' @keywords scRNA-seq
#' @export
#' @examples
#' ComputeModuleEigengene(pbmc)
ComputeModuleEigengene <- function(
  seurat_obj,
  cur_mod,
  modules,
  group.by.vars=NULL,
  verbose=TRUE,
  vars.to.regress = NULL,
  scale.model.use = 'linear',
  pc_dim = 1,
  assay = NULL,
  wgcna_name=NULL, ...
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get the assay
  if(is.null(assay)){assay <- DefaultAssay(seurat_obj)}
  if(dim(seurat_obj@assays[[assay]]@data)[1] == 0){
    stop(paste0("Normalized data slot not found in selected assay ", assay))
  }

  # get genes in this module:
  cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name %>% as.character()

  # subset seurat object by these genes only:
  X_dat <- GetAssayData(seurat_obj, slot='data', assay = assay)[cur_genes,]
  if(dim(seurat_obj@assays[[assay]]@counts)[1] == 0){
    X <- X_dat
  } else{
    X <- GetAssayData(seurat_obj, slot='counts', assay = assay)[cur_genes,]
  }

  # create seurat obj with just these genes
  cur_seurat <- CreateSeuratObject(X, assay = assay, meta.data = seurat_obj@meta.data)
  cur_seurat <- SetAssayData(cur_seurat, slot='data', new.data=X_dat, assay=assay)

  # scale the subsetted expression dataset:
  if(is.null(vars.to.regress)){
    cur_seurat <- ScaleData(cur_seurat, features=rownames(cur_seurat), model.use=scale.model.use)
  } else if(all(vars.to.regress %in% colnames(seurat_obj@meta.data))){
    cur_seurat <- ScaleData(cur_seurat, features=rownames(cur_seurat), model.use=scale.model.use, vars.to.regress=vars.to.regress)
  } else{
    stop(paste0("Some variables specified in vars.to.regress are not found in seurat_obj@meta.data"))
  }

  # compute average expression of each gene
  cur_expr <- GetAssayData(cur_seurat, slot='data')
  expr <- Matrix::t(cur_expr)
  averExpr <- Matrix::rowSums(expr) / ncol(expr)

  # run PCA with Seurat function
  cur_pca <- Seurat::RunPCA(
    cur_seurat,
    features = cur_genes,
    reduction.key=paste0('pca', cur_mod),
    verbose=verbose, ...
  )@reductions$pca
  pc <- cur_pca@cell.embeddings[,pc_dim]
  pc_loadings <- cur_pca@feature.loadings[,pc_dim]

  # correlate average expression with eigengene
  pca_cor <- cor(averExpr, pc)

  # run harmony
  if(!is.null(group.by.vars)){

    # add this PCA as its own reduction in the seurat object
    seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
      embeddings = cur_pca@cell.embeddings,
      assay = assay
    )

    cur_harmony <- harmony::RunHarmony(
      seurat_obj,
      group.by.vars=group.by.vars,
      reduction="ME", verbose=verbose, assay.use=assay, ...
    )@reductions$harmony
    ha <- cur_harmony@cell.embeddings[,pc_dim]
    ha_loadings <- cur_pca@feature.loadings[,pc_dim]

    if(pca_cor < 0){
      cur_harmony@cell.embeddings[,pc_dim] <- -ha
      ha_loadings <- -ha_loadings
    }

    # add harmonized PCA as its own reduction in the seurat object
    seurat_obj@reductions$ME_harmony <- Seurat::CreateDimReducObject(
      embeddings = cur_harmony@cell.embeddings,
      assay = assay
    )

    seurat_obj <- SetMELoadings(
      seurat_obj,
      loadings=ha_loadings,
      harmonized=TRUE,
      wgcna_name=wgcna_name
    )

  }

  if(pca_cor < 0){
    cur_pca@cell.embeddings[,pc_dim] <- -pc
    pc_loadings <- -pc_loadings
  }

  # add this PCA as its own reduction in the seurat object
  seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(
    embeddings = cur_pca@cell.embeddings,
    assay = assay
  )

  seurat_obj <- SetMELoadings(
    seurat_obj,
    loadings=pc_loadings,
    harmonized=FALSE,
    wgcna_name=wgcna_name
  )

  # return seurat object
  seurat_obj

}

#' ModuleEigengenes
#'
#' Computes module eigengenes for all WGCNA co-expression modules
#'
#' @param seurat_obj A Seurat object
#' @param modules table containing module / gene assignments, as in GetModules(seurat_obj).
#' @param group.by.vars groups to harmonize by
#' @param verbose logical indicating whether to print messages
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are "linear", "poisson", or "negbinom"
#' @param pc_dim Which PC to use as the module eigengene? Default to 1.
#' @param assay Assay in seurat_obj to compute module eigengenes. Default is DefaultAssay(seurat_obj)
#' @param exclude_grey logical determining whether to compute MEs for the grey module
#' @param wgcna_name name of the WGCNA experiment
#' @keywords scRNA-seq
#' @export
#' @examples
#'  ModuleEigengenes(pbmc)
ModuleEigengenes <- function(
  seurat_obj,
  group.by.vars=NULL,
  modules=NULL,
  vars.to.regress = NULL,
  scale.model.use = 'linear',
  verbose=TRUE,
  assay = NULL,
  pc_dim = 1,
  exclude_grey = FALSE,
  wgcna_name=NULL, ...
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # are we going to run Harmony?
  harmonized = !is.null(group.by.vars)

  if(harmonized & !any(grepl("ScaleData", seurat_obj@commands))){
      stop('Need to run ScaleData before running ModuleEigengenes with group.by.vars option.')
  }

  # get the assay
  if(is.null(assay)){assay <- DefaultAssay(seurat_obj)}

  # check for the data slot in this assay
  if(dim(seurat_obj@assays[[assay]]@data)[1] == 0){
    stop(paste0("Normalized data slot not found in selected assay ", assay))
  }

  # exclude grey doesn't work yet:
  if(exclude_grey){exclude_grey <- FALSE}


  me_list <- list()
  harmonized_me_list <- list()

  # re-set feature loadings:
  seurat_obj <- SetMELoadings(
    seurat_obj,
    loadings=c(""),
    harmonized=FALSE,
    wgcna_name=wgcna_name
  )
  if(harmonized){
    seurat_obj <- SetMELoadings(
      seurat_obj,
      loadings=c(""),
      harmonized=TRUE,
      wgcna_name=wgcna_name
    )
  }

  # get modules from Seurat object, else use provided modules
  if(is.null(modules)){
    modules <- GetModules(seurat_obj, wgcna_name)
    projected <- FALSE
  } else{
    projected <- TRUE
  }

  # get list of modules:
  mods <- levels(modules$module)
  mods_loop <- mods

  # exclude grey:
  if(exclude_grey){
    mods_loop <- mods[mods != 'grey']
  }

  # loop over modules:
  for(cur_mod in mods_loop){

    print(cur_mod)

    # compute module eigengenes for this module
    seurat_obj <- ComputeModuleEigengene(
      seurat_obj = seurat_obj,
      cur_mod = cur_mod,
      modules=modules,
      group.by.vars=group.by.vars,
      vars.to.regress = vars.to.regress,
      scale.model.use = scale.model.use,
      verbose=verbose,
      pc_dim = pc_dim,
      assay = assay,
      wgcna_name=wgcna_name,
      ...
    )

    # add module eigengene to ongoing list
    cur_me <- seurat_obj@reductions$ME@cell.embeddings[,pc_dim]
    me_list[[cur_mod]] <- cur_me

    # run harmony
    if(harmonized){
      # add module eigengene to ongoing list
      cur_harmonized_me <- seurat_obj@reductions$ME_harmony@cell.embeddings[,pc_dim]
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
  modules$module <- factor(
    as.character(modules$module),
    levels=mods
  )
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

  # remove temp dim reductions by setting to NULL
  seurat_obj@reductions$ME <- NULL
  seurat_obj@reductions$ME_harmony <- NULL

  # return seurat object
  seurat_obj
}

#' ModuleExprScore
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
#' Computes eigengene-based connectivity (kME) based on module eigengenes and the
#' gene expression in selected cell populations.
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing grouping info, ie clusters or celltypes
#' @param group_name name of the group(s) in group.by to use for kME calculation
#' @param corFnc character string specifying the function to be used to calculate co-expression similarity. Defaults to bicor. Any function returning values
#' @param corOptions 	character string specifying additional arguments to be passed to the function given by corFnc. Use "use = 'p', method = 'spearman'" to obtain Spearman correlation. Use "use = 'p'" to obtain Pearson correlation.
#' @param harmonized logical determining whether to use harmonized MEs for kME calculation
#' @param assay Assay in seurat_obj containing expression information.
#' @param slot Slot in specified, default to normalized 'data' slot.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
ModuleConnectivity <- function(
  seurat_obj,
  group.by = NULL,
  group_name = NULL,
  corFnc = 'bicor',
  corOptions = "use='p'",
  harmonized=TRUE,
  assay = NULL,
  slot = 'data',
  wgcna_name = NULL,
  ...
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get module df, wgcna genes, and wgcna params:
  modules <- GetModules(seurat_obj, wgcna_name)
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  # exclude the grey module:
  #modules <- subset(modules, module != 'grey')
  #MEs <- MEs[,colnames(MEs) != 'grey']

  if(is.null(assay)){assay <- DefaultAssay(seurat_obj)}

  if(!is.null(group.by)){
    cells.use <- seurat_obj@meta.data %>% subset(get(group.by) %in% group_name) %>% rownames
    MEs <- MEs[cells.use,]
  } else{
    cells.use <- colnames(seurat_obj)
  }

  # get expression data:
  exp_mat <- GetAssayData(
    seurat_obj,
    assay=assay,
    slot=slot
  )[genes_use,cells.use]

  datExpr <- t(as.matrix(exp_mat))

  kMEs <- WGCNA::signedKME(
    datExpr,
    MEs,
    outputColumnName = "kME",
    corFnc = corFnc,
    corOptions = corOptions,
    ...
  )

  # add module color to the kMEs table
  modules <- modules[,1:3]
  mods <- levels(modules$module)
  colnames(kMEs) <- colnames(MEs)
  kMEs <- kMEs[,mods]
  colnames(kMEs) <- paste0("kME_", colnames(kMEs))
  kMEs <- cbind(modules, kMEs)

  # update the modules table in the Seurat object:
  seurat_obj <- SetModules(seurat_obj, kMEs, wgcna_name)

  seurat_obj

}

#' RunEnrichr
#'
#' Run Enrichr gene set enrichment tests on hdWGCNA modules
#'
#' @param seurat_obj A Seurat object
#' @param dbs character vector of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
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
#' Performs Fisher's Exact Test for overlap between DEGs and hdWGCNA modules.
#'
#' @param seurat_obj A Seurat object
#' @param deg_df DEG table formatted like the output from Seurat's FindMarkers
#' @param fc_cutoff log fold change cutoff for DEGs to be included in the overlap test
#' @param group_col the name of the column in deg_df containing the cell grouping information
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' OverlapModulesDEGs
OverlapModulesDEGs <- function(
  seurat_obj,
  deg_df,
  fc_cutoff = 0.5,
  group_col = 'cluster',
  wgcna_name = NULL,
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  deg_df$group <- deg_df[,group_col]
  cell_groups <- deg_df$group %>% unique

  # get modules,
  modules <- GetModules(seurat_obj)
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
      c(cur_overlap@odds.ratio, cur_overlap@pval, cur_overlap@Jaccard, length(cur_overlap@intersection))
    })) %>% as.data.frame
    colnames(cur_overlap_df) <- c('odds_ratio', 'pval', 'Jaccard', 'size_intersection')
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
  overlap_df <- overlap_df %>% dplyr::select(c(module, group, color, odds_ratio, pval, fdr, Significance, Jaccard, size_intersection))

  overlap_df
}


#' ProjectModules
#'
#' Project a set of co-expression modules from a reference to a query dataset
#'
#' @param seurat_obj A Seurat object where the modules will be projected to
#' @param seurat_ref A Seurat object containing the co-expression modules to be projected
#' @param modules optionally provide a dataframe containing gene-module information to bypass seurat_ref
#' @param group.by.vars groups to harmonize by
#' @param gene_mapping dataframe to map gene names between genomes
#' @param genome1_col column in gene_mapping containing the genes for seurat_ref (reference)
#' @param genome2_col column in gene_mapping containing the genes for seurat_obj (query)
#' @param overlap_proportion the proportion of genes that must be present in seurat_obj for a given module to be projected. Default = 0.5 (50% of genes)
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are "linear", "poisson", or "negbinom"
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_ref@misc slot
#' @param wgcna_name_proj The name of the hdWGCNA experiment to be created for the projected modules in seurat_obj
#' @keywords scRNA-seq
#' @export
#' @examples
#' ProjectModules
ProjectModules <- function(
  seurat_obj,
  seurat_ref,
  modules=NULL,
  group.by.vars=NULL,
  gene_mapping=NULL, # table mapping genes from species 1 to species 2
  genome1_col=NULL, genome2_col=NULL,
  overlap_proportion = 0.5,
  vars.to.regress = NULL,
  scale.model.use = 'linear',
  wgcna_name=NULL, wgcna_name_proj=NULL,
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_ref@misc$active_wgcna}

  # get modules to be projected:
  if(is.null(modules)){
    modules <- GetModules(seurat_ref, wgcna_name)
  } else{
    if(!all(c("gene_name", "module", "color") %in% colnames(modules))){
      stop('Missing columns in modules table. Required columns are gene_name, module, color')
    }
  }

  # cast "modules" to a factor
  if(!is.factor(modules$module)){
    modules$module <- as.factor(modules$module)
  }

  # are we mapping gene names?
  if(!is.null(gene_mapping)){
    modules <- TransferModuleGenome(modules, gene_mapping, genome1_col, genome2_col)
  }

  # what is the proportion of genes in each module that are in seurat_obj?
  mods <- as.character(unique(modules$module))
  mod_props <- unlist(lapply(mods, function(x){
    cur_mod <- subset(modules, module == x)
    sum(cur_mod$gene_name %in% rownames(seurat_obj)) / nrow(cur_mod)
  }))
  mods_keep <- mods[mod_props >= overlap_proportion]
  mods_remove <- mods[mod_props < overlap_proportion]

  if(length(mods) > 0){
    warning(paste0("The following modules will not be projected because too few genes are present in seurat_obj: ", paste(mods_remove, collapse=', ')))
  }

  # only keep modules that have enough overlapping genes
  modules <- subset(modules, module %in% mods_keep) %>%
    dplyr::mutate(module = droplevels(module))

  # get genes that overlap between WGCNA genes & seurat_obj genes:
  gene_names <- modules$gene_name
  genes_use <- intersect(gene_names, rownames(seurat_obj))

  print('n genes:')
  print(length(genes_use))
  print(head(gene_names))

  # subset modules by genes in this seurat object:
  modules <- modules %>% subset(gene_name %in% genes_use)
  print(head(modules))

  # setup new seurat obj for wgcna:
  if(!(wgcna_name_proj %in% names(seurat_obj@misc))){
    seurat_obj <- SetupForWGCNA(
      seurat_obj,
      wgcna_name = wgcna_name_proj,
      features = genes_use
    )
  }

  # project modules:
  print('project modules')
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars=group.by.vars,
    modules = modules,
    vars.to.regress = vars.to.regress,
    scale.model.use = scale.model.use,
    wgcna_name = wgcna_name_proj,
    ...
  )
  print('done')

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

  # TODO
  # need to add a check that there's a one to one mapping
  if(!all(table(gene_mapping[[genome2_col]]) == 1)){
    stop("There can only be one value for each feature in genome2_col in the gene_mapping table.")
  }

  # switch gene names to human:
  gene_mapping <- gene_mapping[,c(genome1_col, genome2_col)]

  # only keep genome1 genes that are in the WGCNA gene list:
  gene_list <- modules$gene_name
  gene_mapping <- gene_mapping[gene_mapping[,genome1_col] %in% gene_list,]

  # subset modules
  modules <- subset(modules, gene_name %in% gene_mapping[,genome1_col])

  # match the order of the given gene list, and remove NA entries
  gene_match <- match(gene_list, gene_mapping[,genome1_col])
  gene_mapping <- na.omit(gene_mapping[gene_match,])

  # update modules table with the new gene names
  # modules <- na.omit(modules[gene_match,])

  print(head(gene_mapping))

  modules$gene_name <- gene_mapping[,genome2_col]
  rownames(modules) <- modules$gene_name

  modules

}


#' ComputeROC
#'
#'
#' @keywords scRNA-seq
#' @export
#' @examples
#' ComputeROC
ComputeROC <- function(
  seurat_obj,
  group.by=NULL, # this needs to be a factor!!!
  split_col=NULL, # needs to be a logical!!!
  features = 'hMEs',
  seurat_test=NULL,
  harmony_group_vars=NULL,
  scale_genes=TRUE,
  verbose=FALSE,
  exp_thresh = 0.75,
  return_seurat=TRUE,
  wgcna_name=NULL, wgcna_name_test=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get Modules
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # if group.by column is null, use Idents:
  if(is.null(group.by)){
    group.by <- 'roc_group'
    seurat_obj@meta.data[[group.by]] <- Idents(seurat_obj)
  }

  # get names of different cell groupings
  groups <- levels(seurat_obj@meta.data[[group.by]])
  groups <- groups[order(groups)]

  # split seurat object into training & testing if the testing is not provided
  if(is.null(seurat_test)){

    print('splitting seurat obj')

    # split into two seurat objects based on train & test split:
    seurat_train <- seurat_obj[,seurat_obj@meta.data[[split_col]]]
    seurat_test <- seurat_obj[,!seurat_obj@meta.data[[split_col]]]

    # project modules for train
    wgcna_name_train <- "ROC"
    seurat_train <- ProjectModules(
      seurat_train,
      seurat_ref = seurat_obj,
      group.by.vars=harmony_group_vars,
      wgcna_name_proj=wgcna_name_train,
      scale_genes=scale_genes, verbose=verbose
    )

  } else{

    if(is.null(wgcna_name_test)){wgcna_name_test <- seurat_test@misc$active_wgcna}

    # if group.by column is null, use Idents:
    if(is.null(group.by)){
      group.by <- 'roc_group'
      seurat_test@meta.data[[group.by]] <- Idents(seurat_test)
    }

    seurat_train <- seurat_obj
    wgcna_name_train <- wgcna_name

  }

  # get names of different cell groupings in test
  groups_test <- levels(seurat_test@meta.data[[group.by]])
  groups_test <- groups_test[order(groups_test)]

  # groups that are in common (some might not have been predicted):
  groups_common <- intersect(
    groups,
    as.character(unique(seurat_test@meta.data[[group.by]]))
  )

  # check if groups are equal:
  if(sum(groups == groups_test) != length(groups)){
    stop("Different groups present in train & test data. Idents likely do not match.")
  }

  # project modules for test if they haven't already been projected:
  if(is.null(GetModules(seurat_test, wgcna_name=wgcna_name_test))){
    wgcna_name_test <- "ROC"
    seurat_test <- ProjectModules(
      seurat_test,
      seurat_ref = seurat_obj,
      group.by.vars=harmony_group_vars,
      wgcna_name_proj=wgcna_name_test,
      scale_genes=scale_genes, verbose=verbose
    )
  }


  # get MEs from seurat object
  if(features == 'hMEs'){
    MEs <- GetMEs(seurat_train, TRUE, wgcna_name_train)
    MEs_p <- GetMEs(seurat_test, TRUE, wgcna_name_test)
  } else if(features == 'MEs'){
    MEs <- GetMEs(seurat_train, FALSE, wgcna_name_train)
    MEs_p <- GetMEs(seurat_test, FALSE, wgcna_name_test)
  } else if(features == 'scores'){
    MEs <- GetModuleScores(seurat_train, wgcna_name_train)
    MEs_p <- GetModuleScores(seurat_test, wgcna_name_test)
    stop("Haven't implemented this one yet >.<")
  } else(
    stop('Invalid feature selection. Valid choices: hMEs, MEs, scores.')
  )


  # add group column to MEs:
  MEs <- as.data.frame(MEs)
  MEs$group <- seurat_train@meta.data[[group.by]]
  MEs_p <- as.data.frame(MEs_p)
  MEs_p$group <- seurat_test@meta.data[[group.by]]


  # only keep common groups:
  MEs <- subset(MEs, group %in% groups_common)
  MEs_p <- subset(MEs_p, group %in% groups_common)

  # compute average MEs in each group:
  avg_MEs <- MEs %>% group_by(group) %>% summarise(across(!!mods, mean))
  avg_MEs_p <- MEs_p %>% group_by(group) %>% summarise(across(!!mods, mean))
  groups <- avg_MEs$group

  # scale each column between 0 & 1
  avg_MEs <- avg_MEs %>% summarise(across(!!mods, scale01))
  avg_MEs_p <- avg_MEs_p %>% summarise(across(!!mods, scale01))

  # convert to binary labels:
  labels <- avg_MEs %>% purrr::map(~ifelse(. >= exp_thresh, TRUE, FALSE))
  labels <- as.data.frame(do.call(cbind, labels));
  rownames(labels) <- as.character(groups)

  # loop over modules to compute ROC curves:
  plot_df <- data.frame()
  conf_df <- data.frame()
  auc_list <- list()
  mod_colors <- list()
  roc_list <- list()
  for(cur_mod in mods){
    print(cur_mod)
    cur_color <- modules %>% subset(module == cur_mod) %>% .$color %>% unique
    mod_colors[[cur_mod]] <- cur_color

    # compute ROC
    rocobj <- pROC::roc(labels[,cur_mod], avg_MEs_p[[cur_mod]])
    auc_list[[cur_mod]] <- as.numeric(rocobj$auc)
    roc_list[[cur_mod]] <- rocobj

    # update plotting df
    cur_df <- data.frame(
      specificity = 1-rocobj$specificities,
      sensitivity = rocobj$sensitivities,
      module = cur_mod,
      color = cur_color,
      auc = as.numeric(rocobj$auc)
    )
    plot_df <- rbind(plot_df, cur_df)

    # update confidence interval df
    cur_conf <- as.data.frame(pROC::ci.se(rocobj))
    cur_conf$sensitivity <- 1-as.numeric(rownames(cur_conf))
    cur_conf$module <- cur_mod
    cur_conf$color <- cur_color
    conf_df <- rbind(conf_df, cur_conf)

  }
  colnames(conf_df)[1:3] <- c('lo', 'mid', 'hi')

  # return ROC tables & objects:
  roc_info = list(
    roc = plot_df,
    conf = conf_df,
    objects = roc_list
  )

  if(return_seurat){
    # add roc info to seurat object:
    seurat_obj <- SetROCData(seurat_obj, roc_info, wgcna_name)
    out <- seurat_obj
  } else{
    out <- roc_info
  }

  out
}


#' Scan gene promoters for a set of TF PWMs
#'
#'
#' @keywords scRNA-seq
#' @export
#' @examples
#' MotifScan
MotifScan <- function(
  seurat_obj,
  species_genome, # hg38, mm10, etc...
  pfm, # matrix set from JASPAR2020 for example
  EnsDb, # Ensembl database such as EnsDb.Mmusculus.v79
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get a dataframe of just the motif name and the motif ID:
  motif_df <- data.frame(
    motif_name = purrr::map(1:length(pfm), function(i){pfm[[i]]@name}) %>% unlist,
    motif_ID = purrr::map(1:length(pfm), function(i){pfm[[i]]@ID}) %>% unlist
  )

  # get promoters of protein coding genes from the given Ensdb
  # note: everything breaks if I try to use X & Y chromosomes.
  gene.promoters <- ensembldb::promoters(EnsDb, filter = ~ gene_biotype == "protein_coding") %>%
    subset(seqnames %in% c(1:100))
  gene.coords <- ensembldb::genes(EnsDb, filter = ~ gene_biotype == "protein_coding") %>%
    subset(seqnames %in% c(1:100))

  # add the gene name to the promoter object
  gene.promoters$symbol <- gene.coords$symbol[match(gene.promoters$gene_id, names(gene.coords))]

  # drop unnecessary chromosomes
  gene.promoters <- keepSeqlevels(gene.promoters, value= levels(droplevels(seqnames(gene.promoters))))

  # rename seqlevels to add 'chr', remove X&Y chromosomes because they break everything
  old_levels <- levels(seqnames(gene.promoters))
  new_levels <- ifelse(old_levels %in% c('X', 'Y'), old_levels, paste0('chr', old_levels))
  gene.promoters <- renameSeqlevels(gene.promoters, new_levels)

  # set the genome (not sure if we NEED to do this...)
  genome(seqinfo(gene.promoters)) <- species_genome

  # set up promoters object that only has the necessary info for motifmatchr
  my_promoters <- GRanges(
    seqnames =  droplevels(seqnames(gene.promoters)),
    IRanges(
      start = start(gene.promoters),
      end = end(gene.promoters)
    ),
    symbol = gene.promoters$symbol,
    genome=species_genome
  )

  # scan these promoters for motifs:
  print('Matching motifs...')
  motif_ix <- motifmatchr::matchMotifs(pfm, my_promoters, genome=species_genome)

  # get the matches
  tf_match <- motifmatchr::motifMatches(motif_ix)
  rownames(tf_match) <- my_promoters$symbol

  # use motif names as the column names:
  colnames(tf_match) <- motif_df$motif_name

  # only keep genes that are in the Seurat object and in the given EnsDb:
  gene_list <- rownames(seurat_obj)
  gene_list <- gene_list[gene_list %in% rownames(tf_match)]
  tf_match <- tf_match[gene_list,]

  # get list of target genes for each TF:
  print('Getting TF target genes...')
  tfs <- motif_df$motif_name
  tf_targets <- list()
  n_targets <- list()
  for(cur_tf in tfs){
    tf_targets[[cur_tf]] <- names(tf_match[,cur_tf][tf_match[,cur_tf]])
    n_targets[[cur_tf]] <- length(tf_targets[[cur_tf]] )
  }
  n_targets <- unlist(n_targets)

  # add number of target genes to motif df
  motif_df$n_targets <- n_targets

  # add info to seurat object
  seurat_obj <- SetMotifMatrix(seurat_obj, tf_match)
  seurat_obj <- SetMotifs(seurat_obj, motif_df)
  seurat_obj <- SetMotifTargets(seurat_obj, tf_targets)
  seurat_obj <- SetPFMList(seurat_obj, pfm)

  seurat_obj
}


#' Overlap modules with TF target genes
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' OverlapModulesDEGs
OverlapModulesMotifs <- function(
  seurat_obj, wgcna_name = NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get tf info
  tf_match <- GetMotifMatrix(seurat_obj)
  tf_targets <- GetMotifTargets(seurat_obj)
  motif_df <- GetMotifs(seurat_obj)

  # size of genome based on # genes in Seurat object:
  genome.size <- nrow(seurat_obj)

  cur_mod <- mods[1]
  cur_module_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name
  cur_targets <- as.character(unlist(tf_targets[1]))

  overlap_df <- do.call(rbind, lapply(mods, function(cur_mod){
    cur_module_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name
    cur_overlap_df <- do.call(rbind, lapply(tf_targets, function(cur_targets){
      cur_overlap <- testGeneOverlap(newGeneOverlap(
          cur_module_genes,
          as.character(unlist(cur_targets)),
          genome.size=genome.size
      ))
      c(cur_overlap@odds.ratio, cur_overlap@pval, cur_overlap@Jaccard, length(cur_overlap@intersection))
    })) %>% as.data.frame

    colnames(cur_overlap_df) <- c('odds_ratio', 'pval', 'Jaccard', 'size_intersection')
    cur_overlap_df$module <- cur_mod
    cur_overlap_df$tf <- names(tf_targets)

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
  overlap_df <- overlap_df %>% dplyr::select(c(module, tf, color, odds_ratio, pval, fdr, Significance, Jaccard, size_intersection))

  # add overlap df to Seurat obj
  seurat_obj <- SetMotifOverlap(seurat_obj, overlap_df, wgcna_name)

  seurat_obj

}


#' MotifTargetScore
#'
#' Computes gene expression scores for TF Motif target genes based on the MotifScan.
#'
#' @param seurat_obj A Seurat object
#' @param method Seurat or UCell?
#' @keywords scRNA-seq
#' @export
#' @examples
#' MotifTargetScore(pbmc)
MotifTargetScore <- function(
  seurat_obj,
  method='Seurat',
  wgcna_genes=TRUE,
  wgcna_name=NULL,
  ...
 ){

   if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

   # get modules:
   modules <- GetModules(seurat_obj, wgcna_name)

   # get TF target genes:
   target_genes <- GetMotifTargets(seurat_obj)

   # subset by WGCNA genes only:
   if(wgcna_genes){
     target_genes <- lapply(target_genes, function(x){
       x[x %in% modules$gene_name]
     })
   }

   # run gene scoring function
   if(method == "Seurat"){
     tf_scores <- Seurat::AddModuleScore(
       seurat_obj, features=target_genes, ...
     )@meta.data
   } else if(method == "UCell"){
     tf_scores <- UCell::AddModuleScore_UCell(
       seurat_obj, features=target_genes, ...
     )@meta.data
   } else{
     stop("Invalid method selection. Valid choices are Seurat, UCell")
   }

   tf_scores <- tf_scores[,(ncol(tf_scores)-length(target_genes)+1):ncol(tf_scores)]
   colnames(tf_scores) <- names(target_genes)

   # add tf scores to seurat object:
   seurat_obj <- SetMotifScores(seurat_obj, tf_scores, wgcna_name)

   seurat_obj

 }


#' RunModuleUMAP
#'
#' Run UMAP on co-expression matrix using hub genes as features.
#'
#' @param seurat_obj A Seurat object
#' @param n_hubs number of hub genes to use in the UMAP computation
#' @param exclude_grey logical indicating whether to include grey module
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Additional parameters supplied to uwot::umap
#' @keywords scRNA-seq
#' @export
#' @examples
RunModuleUMAP <- function(
  seurat_obj,
  features = "TOM", # "TOM" or "kME"
  n_hubs = 10,
  exclude_grey = TRUE,
  wgcna_name = NULL,
  n_neighbors= 25,
  metric = "cosine",
  spread=1,
  min_dist=0.4,
  supervised = FALSE,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules,
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)

  # check if we have eigengene-based connectivities:
  if(!all(paste0('kME_', as.character(mods)) %in% colnames(modules))){
    stop('Eigengene-based connectivity (kME) not found. Did you run ModuleEigengenes and ModuleConnectivity?')
  }

  if(exclude_grey){
    mods <- mods[mods != 'grey']
    modules <- subset(modules, module != 'grey')
  }

  # get hub genes:
  hub_list <- lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(n_hubs) %>% .$gene_name
  })
  names(hub_list) <- mods

  # get all genes that aren't in gray mod
  selected_genes <- modules[modules$module %in% mods,'gene_name']

  # get the TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # subset the TOM for umap
  # keep all genes as rows, and keep only hubs as cols (features)
  feature_mat <- TOM[selected_genes,unlist(hub_list)]

  # run UMAP
  if(supervised){
    print('running supervised UMAP:')
    hub_umap <-  uwot::umap(
      X = feature_mat,
      min_dist = min_dist,
      n_neighbors= n_neighbors,
      metric = metric,
      spread=spread,
      y = modules$module, # for supervised UMAP
      ...
    )
  } else {
    hub_umap <-  uwot::umap(
      X = feature_mat,
      min_dist = min_dist,
      n_neighbors= n_neighbors,
      metric = metric,
      spread=spread,
      ...
    )
  }


  # set up plotting df
  plot_df <- as.data.frame(hub_umap)
  colnames(plot_df) <- c("UMAP1", "UMAP2")
  plot_df$gene <- rownames(feature_mat)

  # add module color, and hub gene status to the plotting df:
  ix <- match(plot_df$gene, modules$gene_name)
  plot_df$module <- modules$module[ix]
  plot_df$color <- modules$color[ix]
  plot_df$hub <- ifelse(
    plot_df$gene %in% as.character(unlist(hub_list)), 'hub', 'other'
  )

  # get kME values for each gene
  kMEs <- do.call(rbind, lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur <- cur[,c('gene_name', paste0('kME_', cur_mod))]
    colnames(cur) <- c('gene_name', 'kME')

    # scale kMEs between 0 & 1:
    cur$kME <- scale01(cur$kME)
    cur
  }))
  ix <- kMEs$gene_name[match(plot_df$gene, kMEs$gene_name)]
  plot_df$kME <- kMEs[ix, 'kME']

  # add the UMAP to the Seurat object:
  seurat_obj <- SetModuleUMAP(seurat_obj, plot_df, wgcna_name)

  seurat_obj
}


#' ModuleTraitCorrelation'
#'
#' Correlates categorical and numeric variables with Module Eigengenes or hub-gene scores.
#'
#'
#' @param seurat_obj A Seurat object
#' @param seurat_obj A list of column names in the Seurat object's metadata that you wish to correlate with each module.
#' Traits must be a categorical variable (not a character vector), or a numeric variable.
#' @param features Which features to use to summarize each modules? Valid choices are hMEs, MEs, or scores
#' @param cor_meth Which method to use for correlation? Valid choices are pearson, spearman, kendall.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModuleTraitCorrelation
ModuleTraitCorrelation <- function(
  seurat_obj,
  traits,
  group.by = NULL,
  features = 'hMEs',
  cor_method = 'pearson',
  subset_by = NULL,
  subset_groups = NULL,
  wgcna_name = NULL,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get MEs, module data from seurat object
  if(features == 'hMEs'){
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  } else if(features == 'MEs'){
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  } else if(features == 'scores'){
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  } else{
    stop('Invalid feature selection. Valid choices: hMEs, MEs, scores, average')
  }

  # subset?
  if(!is.null(subset_by)){
    print('subsetting')
    seurat_full <- seurat_obj
    MEs <- MEs[seurat_obj@meta.data[[subset_by]] %in% subset_groups,]
    seurat_obj <- seurat_obj[,seurat_obj@meta.data[[subset_by]] %in% subset_groups]
  }

  # check if traits are in the seurat object:
  if(sum(traits %in% colnames(seurat_obj@meta.data)) != length(traits)){
    stop(paste('Some of the provided traits were not found in the Seurat obj:', paste(traits[!(traits %in% colnames(seurat_obj@meta.data))], collapse=', ')))
  }

  #use idents as grouping variable if not specified
  if(is.null(group.by)){
    group.by <- 'temp_ident'
    seurat_obj$temp_ident <- Idents(seurat_obj)
  }

  # check the class of each trait provided:
  valid_types <- c('numeric', 'factor', 'integer')
  data_types <- sapply(traits, function(x){class(seurat_obj@meta.data[,x])})

  if(!all(data_types %in% valid_types)){
    incorrect <- traits[!(data_types %in% valid_types)]
    stop(paste0('Invalid data types for ', paste(incorrect, collapse=', '), '. Accepted data types are numeric, factor, integer.'))
  }

  # print warnings about factor levels:
  if(any(data_types == 'factor')){
    factor_traits <- traits[data_types == 'factor']
    for(tr in factor_traits){
      warning(paste0("Trait ", tr, ' is a factor with levels ', paste0(levels(seurat_obj@meta.data[,tr]), collapse=', '), '. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?'))
    }
  }

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # get trait table:
  trait_df <- seurat_obj@meta.data[,traits]

  # cast vector to data frame if there's only one trait
  if(length(traits == 1)){
    trait_df <- data.frame(x = trait_df)
    colnames(trait_df) <- traits
  }

  # convert factors to numeric
  if(any(data_types == 'factor')){
    factor_traits <- traits[data_types == 'factor']
    for(tr in factor_traits){
      trait_df[,tr] <- as.numeric(trait_df[,tr])
    }
  }

  # correlate all cells:
  cor_list <- list(); pval_list <- list(); fdr_list <- list()

  # correlation:
  temp <- Hmisc::rcorr(as.matrix(trait_df), as.matrix(MEs), type=cor_method)

  # get the coefficient & p-val
  cur_cor <- temp$r[traits,mods]
  cur_p <- temp$P[traits,mods]

  # compute FDR:
  p_df <- cur_p %>%
    reshape2::melt()

    if(length(traits) == 1){

      tmp <- rep(mods, length(traits))
      tmp <- factor(tmp, levels = mods)
      tmp <- tmp[order(tmp)]

      p_df$Var1 <- traits
      p_df$Var2 <- tmp
      rownames(p_df) <- 1:nrow(p_df)
      p_df <- dplyr::select(p_df, c(Var1, Var2, value))
    }

  p_df <- p_df %>%
  dplyr::mutate(fdr=p.adjust(value, method='fdr')) %>%
  dplyr::select(c(Var1, Var2, fdr))

  # reshape to match cor & pval
  cur_fdr <- reshape2::dcast(p_df, Var1 ~ Var2, value.var='fdr')
  rownames(cur_fdr) <- cur_fdr$Var1
  cur_fdr <- cur_fdr[,-1]

  # add to list
  cor_list[["all_cells"]] <- cur_cor
  pval_list[["all_cells"]] <- cur_p
  fdr_list[["all_cells"]] <- cur_fdr

  trait_df <- cbind(trait_df, seurat_obj@meta.data[,group.by])
  colnames(trait_df)[ncol(trait_df)] <- 'group'

  MEs <- cbind(as.data.frame(MEs), seurat_obj@meta.data[,group.by])
  colnames(MEs)[ncol(MEs)] <- 'group'

  if(class(seurat_obj@meta.data[,group.by]) == 'factor'){
    group_names <- levels(seurat_obj@meta.data[,group.by])
  } else{
    group_names <- levels(as.factor(seurat_obj@meta.data[,group.by]))
  }

  # do the correlation for each group:
  trait_list <- dplyr::group_split(trait_df, group, .keep=FALSE)
  ME_list <- dplyr::group_split(MEs, group, .keep=FALSE)
  names(trait_list) <- group_names
  names(ME_list) <- group_names

  for(i in names(trait_list)){
    # cor_list[[i]] <- cor(as.matrix(trait_list[[i]]), as.matrix(ME_list[[i]]), method=cor_method)

    # testing other correlation function:
    temp <- Hmisc::rcorr(as.matrix(trait_list[[i]]), as.matrix(ME_list[[i]]))

    cur_cor <- temp$r[traits,mods]
    cur_p <- temp$P[traits,mods]

    # compute FDR:
    p_df <- cur_p %>%
      reshape2::melt()


    if(length(traits) == 1){

      tmp <- rep(mods, length(traits))
      tmp <- factor(tmp, levels = mods)
      tmp <- tmp[order(tmp)]

      p_df$Var1 <- traits
      p_df$Var2 <- tmp
      rownames(p_df) <- 1:nrow(p_df)
      p_df <- dplyr::select(p_df, c(Var1, Var2, value))
    }

    p_df <- p_df %>%
      dplyr::mutate(fdr=p.adjust(value, method='fdr')) %>%
      dplyr::select(c(Var1, Var2, fdr))

    # reshape to match cor & pval
    cur_fdr <- reshape2::dcast(p_df, Var1 ~ Var2, value.var='fdr')
    rownames(cur_fdr) <- cur_fdr$Var1
    cur_fdr <- cur_fdr[,-1]

    # add to list
    cor_list[[i]] <- cur_cor
    pval_list[[i]] <- cur_p
    fdr_list[[i]] <- as.matrix(cur_fdr)

  }

  # add Module-trait correlations to the seruat object:
  mt_cor <- list(
    'cor' = cor_list,
    'pval' = pval_list,
    'fdr' = fdr_list
  )

  if(!is.null(subset_by)){
    seurat_full <- SetModuleTraitCorrelation(seurat_full, mt_cor, wgcna_name)
    seurat_obj <- seurat_full
  } else{
    seurat_obj <- SetModuleTraitCorrelation(seurat_obj, mt_cor, wgcna_name)
  }

  seurat_obj
}




#' ModulePreservation
#'
#' Computes module preservation statistics in a query dataset for a given reference dataset
#'
#' @param seurat_obj A Seurat object
#' @param seurat_ref A Seurat object serving as the reference for the module preservation analysis
#' @param name The name to give the module preservation analysis.
#' @param n_permutations Number of permutations for the module preservation test.
#' @param parallel logical determining whether to run preservation analysis in parallel
#' @param seed random seed for the permutation analysis.
#' @param return_raw if TRUE, returns the module preservation statistics, else returns seurat_obj with the stats added to the hdWGCNA experiment.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param wgcna_name_ref The name of the hdWGCNA experiment in the seurat_ref@misc slot
#' @keywords scRNA-seq
#' @export
#' @examples
#' ModulePreservation
ModulePreservation <- function(
  seurat_obj,
  seurat_ref,
  name,
  n_permutations = 500,
  parallel = FALSE,
  seed = 12345,
  gene_mapping = NULL,
  genome1_col = NULL,
  genome2_col = NULL,
  return_raw = FALSE,
  wgcna_name = NULL,
  wgcna_name_ref = NULL,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  if(is.null(wgcna_name_ref)){wgcna_name_ref <- seurat_ref@misc$active_wgcna}

  print('Setup datasets')

  # get datExpr for reference and query:
  datExpr_ref <- GetDatExpr(seurat_ref, wgcna_name_ref)
  datExpr_query <- GetDatExpr(seurat_obj, wgcna_name)

  # change the gene names to match:
  if(!is.null(gene_mapping)){
    gene_match <- match(colnames(datExpr_query), gene_mapping[,genome2_col])
    gene_mapping <- na.omit(gene_mapping[gene_match,])
    colnames(datExpr_query)  <- gene_mapping[,genome1_col]

    print(head(colnames(datExpr_query)))
    print(head(GetModules(seurat_obj)$gene_name))
  }

  # set up multiExpr:
  setLabels <- c("ref", "query")
  multiExpr <- list(
    ref = list(data=datExpr_ref),
    query = list(data=datExpr_query)
  )
  ref_modules <- list(ref = GetModules(seurat_ref)$module)

  print('ref:')
  print(dim(datExpr_ref))
  print('query:')
  print(dim(datExpr_query))

  print('Run Module Preservation')

  # run the module preservation test:
  mp <- WGCNA::modulePreservation(
    multiExpr,
    ref_modules,
    referenceNetworks = 1,
    nPermutations = n_permutations,
    randomSeed = seed,
    quickCor = 0,
    parallelCalculation = parallel,
    ...
  )

  if(return_raw){return(mp)}

  # get the stats obs and stats Z tables
  ref <- 1; test <- 2;
  statsObs <- cbind(
    mp$quality$observed[[ref]][[test]],
    mp$preservation$observed[[ref]][[test]]
  )
  statsZ <- cbind(
    mp$quality$Z[[ref]][[test]],
    mp$preservation$Z[[ref]][[test]]
  )

  # add stats to the seurat object:
  mod_pres <- list('obs'  = statsObs, 'Z' = statsZ)
  seurat_obj <- SetModulePreservation(seurat_obj, mod_pres, name, wgcna_name)

  seurat_obj

}


#' FindAllDMEs
#'
#' Modification of the Seurat function FindMarkers used to perform iterative one-versus-all differential module eigengene testing given a group of cells (ie clusters, cell types, etc).
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing cell grouping information
#' @param harmonized logical determining whether or not to use harmonized MEs
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param add_missing logical determining whether or not to add missing modules back into the resulting dataframe with NA values.
#' @param ... Additional parameters for the FindMarkers function
#' @keywords scRNA-seq
#' @export
#' @return A dataframe contaning differential ME results
#' @examples
#' FindAllDMEs
FindAllDMEs <- function(
  seurat_obj,
  group.by,
  harmonized=TRUE,
  wgcna_name=NULL,
  add_missing=FALSE,
  test.use='wilcox',
  only.pos=FALSE,
  logfc.threshold = 0,
  min.pct=0,
  verbose=FALSE,
  pseudocount.use=0,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get list of groups
  groups <- as.character(unique(seurat_obj@meta.data[[group.by]]))

  # get module eigengenes, remove grey module
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  MEs <- MEs[, colnames(MEs) != 'grey']

  # set all negative values to zero
  MEs[MEs < 0] <- 0

  # transpose MEs so cells are columns and rows are MEs
  MEs <- t(MEs)

  # create new assay for MEs:
  ME_assay <- Seurat::CreateAssayObject(MEs)

  DMEs_list <- list()
  for(cur_group in groups){

    print(cur_group)

    barcodes1 <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] == cur_group]
    barcodes2 <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] != cur_group]

    # run differential test on the ME assay
    DMEs <- FindMarkers(
      ME_assay,
      cells.1 = barcodes1,
      cells.2 = barcodes2,
      slot='counts',
      test.use=test.use,
      only.pos=only.pos,
      logfc.threshold=logfc.threshold,
      min.pct=min.pct,
      verbose=verbose,
      pseudocount.use=pseudocount.use,
      ...
    )
    DMEs$module <- rownames(DMEs)
    DMEs$group <- cur_group

    # add missing modules if specified
    if(add_missing){
      missing_mods <- rownames(MEs)[!(rownames(MEs) %in% DMEs$module)]
      for(cur_mod in missing_mods){
        DMEs[cur_mod,] <- NA
        DMEs[cur_mod,'module'] <- cur_mod
      }
    }

    # add to ongoing list:
    DMEs_list[[cur_group]] <- DMEs

  }

  DMEs <- do.call(rbind, DMEs_list)
  DMEs
}


#' FindDMEs
#'
#' Modification of the Seurat function FindMarkers used to perform differential module eigengene testing between two sets of cell barcodes.
#'
#' @param seurat_obj A Seurat object
#' @param barcodes1 character vector containing cell barcodes for the first group to test. Positive fold-change means up-regulated in this group.
#' @param barcodes2 character vector containing cell barcodes for the second group to test. Negative fold-change means up-regulated in this group.
#' @param harmonized logical determining whether or not to use harmonized MEs
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param add_missing logical determining whether or not to add missing modules back into the resulting dataframe with NA values.
#' @param ... Additional parameters for the FindMarkers function
#' @keywords scRNA-seq
#' @export
#' @return A dataframe contaning differential ME results
#' @examples
#' FindDMEs
FindDMEs <- function(
  seurat_obj,
  barcodes1,
  barcodes2,
  harmonized=TRUE,
  wgcna_name=NULL,
  add_missing=FALSE,
  test.use='wilcox',
  only.pos=FALSE,
  logfc.threshold = 0,
  min.pct=0,
  verbose=FALSE,
  pseudocount.use=0,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # ensure that selected barcodes are in the seurat obj
  if(!(all(barcodes1 %in% colnames(seurat_obj)))){
    stop('Some barcodes in barcodes1 not found in colnames(seurat_obj).')
  }
  if(!(all(barcodes2 %in% colnames(seurat_obj)))){
    stop('Some barcodes in barcodes2 not found in colnames(seurat_obj).')
  }

  # check for overlap in the two groups of bacrodes:
  if(length(intersect(barcodes1, barcodes2)) > 0){
    stop('Some barcodes overlap in barcodes1 and barcodes2')
  }

  # get module eigengenes, remove grey module
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  MEs <- MEs[, colnames(MEs) != 'grey']
  print(dim(MEs))
  print(colnames(MEs))

  # set all negative values to zero
  MEs[MEs < 0] <- 0

  # transpose MEs so cells are columns and rows are MEs
  MEs <- t(MEs)

  # create new assay for MEs:
  ME_assay <- Seurat::CreateAssayObject(MEs)

  # run differential test on the ME assay
  DMEs <- FindMarkers(
    ME_assay,
    cells.1 = barcodes1,
    cells.2 = barcodes2,
    slot='counts',
    test.use=test.use,
    only.pos=only.pos,
    logfc.threshold=logfc.threshold,
    min.pct=min.pct,
    verbose=verbose,
    pseudocount.use=pseudocount.use,
    ...
  )
  DMEs$module <- rownames(DMEs)

  # add missing modules if specified
  if(add_missing){
    missing_mods <- rownames(MEs)[!(rownames(MEs) %in% DMEs$module)]
    for(cur_mod in missing_mods){
      DMEs[cur_mod,] <- NA
      DMEs[cur_mod,'module'] <- cur_mod
    }
  }

  DMEs

}
