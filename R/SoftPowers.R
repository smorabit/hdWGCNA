
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
