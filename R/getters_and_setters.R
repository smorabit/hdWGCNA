
############################
# Active WGCNA
###########################

#' SetActiveWGCNA
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetActiveWGCNA <- function(seurat_obj, wgcna_name){

  # set the active_wgcna variable
  seurat_obj@misc$active_wgcna <- wgcna_name

  # initialize empty list for this WGCNA if it doesn't exist yet
  if(!(wgcna_name %in% names(seurat_obj@misc))){
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]] <- list()
  }
  seurat_obj
}

#' GetActiveWGCNA
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
GetActiveWGCNA <- function(seurat_obj){
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]
}

#' GetActiveWGCNAName
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
GetActiveWGCNAName <- function(seurat_obj){
  seurat_obj@misc$active_wgcna
}

#' CheckWGCNAName
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
CheckWGCNAName <- function(seurat_obj, wgcna_name){
  wgcna_name %in% names(seurat_obj@misc)
}

# get any WGCNA data, but by default get the active
GetWGCNA <- function(seurat_obj, wgcna_name=NULL){

  # test if wgcna_name is valid (TODO)

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  seurat_obj@misc[[wgcna_name]]
}

############################
# metacell object
###########################

#' SetMetacellObject
#'
#' @param seurat_obj A Seurat object
#' @param metacell_obj metacell Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMetacellObject <- function(seurat_obj, metacell_obj, wgcna_name=NULL){
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # add metacell obj to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj <- metacell_obj
  seurat_obj
}

#' GetMetacellObject
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMetacellObject <- function(seurat_obj,  wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  input_class <- class(seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj)
  if(input_class == "Seurat"){
    return(seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj)
  } else if(input_class == "character") {
    metacell_location <- seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj
    return(seurat_obj@misc[[metacell_location]]$wgcna_metacell_obj)
  } else{
    return(NULL)
  }
}

############################
# WGCNA genes
###########################

#' SetWGCNAGenes
#'
#' @param seurat_obj A Seurat object
#' @param gene_list vector of genes to be used for WGCNA
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetWGCNAGenes <- function(seurat_obj, gene_list, wgcna_name=NULL){

  # test if wgcna_name is valid (TODO)
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # add gene list to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_genes <- gene_list
  seurat_obj
}

#' GetWGCNAGenes
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetWGCNAGenes <- function(seurat_obj, wgcna_name=NULL){
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_genes
}

############################
# datExpr
###########################


#' SetDatExpr
#'
#' This function specifies the gene expression matrix for co-expression network analysis.
#'
#' @param seurat_obj A Seurat object
#' @param group_name A string containing a group present in the provided group.by column or in the Seurat Idents. A character vector can be provided to select multiple groups at a time.
#' @param use_metacells A logical determining if we use the metacells (TRUE) or the full expression matrix (FALSE)
#' @param group.by A string containing the name of a column in the Seurat object with cell groups (clusters, cell types, etc). If NULL (default), hdWGCNA uses the Seurat Idents as the group.
#' @param multi.group.by A string containing the name of a column in the Seurat object with groups for consensus WGCNA (dataset, sample, condition, etc)
#' @param multi_group_name A string containing the name of a group present in the multi.group.by column.
#' @param assay The name of the assay in the Seurat object
#' @param slot The name of the slot in the Seurat object (counts, data)
#' @param mat A Matrix containing gene expression data. Supplying a matrix using this parameter ignores other options. This is almost exclusively used for pseudobulk analysis.
#' @param wgcna_name A string containing the name of the WGCNA slot in seurat_obj@misc. Default = NULL which retrieves the currently active WGCNA data
#' @details
#' SetDatExpr is a critical function of the hdWGCNA pipeline that determines the gene expession 
#' matrix that will be used for network analysis. We typically use this function to select a 
#' cell type or group of cell types for network analysis using the group.by parameter, but we provide 
#' additional parameters for further customization.
#' @keywords scRNA-seq
#' @export
SetDatExpr <- function(
  seurat_obj,
  group_name,
  use_metacells=TRUE,
  group.by=NULL,
  multi.group.by = NULL,
  multi_group_name = NULL,
  return_seurat = TRUE,
  assay=NULL,
  slot = 'data',
  mat=NULL,
  wgcna_name=NULL,
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  if(!CheckWGCNAName(seurat_obj, wgcna_name)){
    stop(paste0("Invalid wgcna_name supplied: ", wgcna_name))
  }  

  if(is.null(assay)){
      assay <- DefaultAssay(s_obj)
      warning(paste0('assay not specified, trying to use assay ', assay))
  }

  # check that selected assay is in the seurat object 
  if(!(assay %in% Assays(seurat_obj))){
    stop(paste0('Invalid choice of assay: ', assay, ' not found in Assays(seurat_obj).'))
  }

  # check that slot is valid 
  if(!(slot %in% c('counts', 'data', 'scale.data'))){
    stop('Invalid choice of slot. Valid choices are counts, data, or scale.data.')
  }

  # check that group.by is valid 
  if(!(group.by %in% names(seurat_obj@meta.data))){
    stop(paste0(group.by, ' not found in seurat_obj@meta.data.'))
  }
  if(!(class(seurat_obj@meta.data[,group.by]) %in% c('character', 'factor'))){
    stop('Selected group.by must be a character or a factor, but ', group.by, ' is a ', class(seurat_obj@meta.data[,group.by]), '.')
  }

    # check that multi.group.by is valid 
  if(!(multi.group.by %in% names(seurat_obj@meta.data))){
    stop(paste0(multi.group.by, ' not found in seurat_obj@meta.data.'))
  }
  if(!(class(seurat_obj@meta.data[,multi.group.by]) %in% c('character', 'factor'))){
    stop('Selected group.by must be a character or a factor, but ', group.by, ' is a ', class(seurat_obj@meta.data[,group.by]), '.')
  }

  # get parameters from seurat object
  params <- GetWGCNAParams(seurat_obj, wgcna_name)
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)

  # was a matrix supplied?
  if(is.null(mat)){
      
    # get metacell object
    m_obj <- GetMetacellObject(seurat_obj, wgcna_name)

    # use metacells or whole seurat object?
    if(use_metacells & !is.null(m_obj)){
      s_obj <- m_obj
    } else{
      if(is.null(m_obj)){warning("Metacell Seurat object not found. Using full Seurat object instead.")}
      s_obj <- seurat_obj
    }

    # get the metadata from the seurat object:
    seurat_meta <- s_obj@meta.data

    # check the group.by params
    if(!is.null(group.by)){

      # check that group.by is in the Seurat object & in the metacell object:
      if(!(group.by %in% colnames(s_obj@meta.data))){
        m_cell_message <- ""
        if(use_metacells){m_cell_message <- "metacell"}
        stop(paste0(group.by, ' not found in the meta data of the ', m_cell_message, ' Seurat object'))
      }

      # check that the selected groups are in the Seurat object:
      if(!all(group_name %in% s_obj@meta.data[[group.by]])){
        groups_not_found <- group_name[!(group_name %in% s_obj@meta.data[[group.by]])]
        stop(
          paste0("Some groups in group_name are not found in the seurat_obj: ", paste(groups_not_found, collapse=', '))
        )
      }
    }

    # check the multi.group.by params
    if(!is.null(multi.group.by)){

      # check that group.by is in the Seurat object & in the metacell object:
      if(!(multi.group.by %in% colnames(s_obj@meta.data))){
        m_cell_message <- ""
        if(use_metacells){m_cell_message <- "metacell"}
        stop(paste0(multi.group.by, ' not found in the meta data of the ', m_cell_message, ' Seurat object'))
      }

      # check that the selected groups are in the Seurat object:
      if(!all(multi_group_name %in% s_obj@meta.data[[multi.group.by]])){
        groups_not_found <- multi_group_name[!(multi_group_name %in% s_obj@meta.data[[multi.group.by]])]
        stop(
          paste0("Some groups in group_name are not found in the seurat_obj: ", paste(groups_not_found, collapse=', '))
        )
      }
    }

    # columns to group by for cluster/celltype
    if(!is.null(group.by)){
      seurat_meta <- seurat_meta %>% subset(get(group.by) %in% group_name)
    }


    # subset further if multiExpr:
    if(!is.null(multi.group.by)){
      seurat_meta <- seurat_meta %>% subset(get(multi.group.by) %in% multi_group_name)
    }

    # get list of cells to use
    cells <- rownames(seurat_meta)

    # get expression data from seurat obj
    datExpr <- as.data.frame(
      Seurat::GetAssayData(
        s_obj,
        assay=assay,
        slot=slot
      )[genes_use,cells]
    )

    # transpose data
    datExpr <- as.data.frame(t(datExpr))
  } else{

    datExpr <- mat

    # cast it to a dataframe
    if(any(class(datExpr) != 'data.frame')){
      datExpr <- as.data.frame(datExpr)
    }

    # are the colnames genes?
    if(!all(colnames(datExpr) %in% rownames(seurat_obj))){
      stop("colnames of the provided matrix are invalid. Make sure that the colnames are features (genes), and that all of these features are in the seurat_obj")
    }

    # subset the datExpr by the WGCNA genes:
    datExpr <- datExpr[,genes_use]

  }

  if(return_seurat){
    gene_list <- genes_use[WGCNA::goodGenes(datExpr, ...)]
    datExpr <- datExpr[,gene_list]

    # update the WGCNA gene list:
    seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)

    # set the datExpr in the Seurat object
    seurat_obj@misc[[wgcna_name]]$datExpr <- datExpr
    out <- seurat_obj
  } else{
    out <- datExpr
  }
  out
}



#' GetDatExpr
#'
#' This function gets the WGCNA expression matrix.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
GetDatExpr <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$datExpr

}


#' SetMultiExpr
#'
#' This function sets up the expression matrix input for consensus WGCNA based on
#' the metacell expression matrix, the full expression matrix, or a provided pseudobulk expression matrix.
#'
#' @param seurat_obj A Seurat object
#' @param group_name A string containing a group present in the provided group.by column or in the Seurat Idents.
#' @param use_metacells A logical determining if we use the metacells (TRUE) or the full expression matrix (FALSE)
#' @param group.by A string containing the name of a column in the Seurat object with cell groups (clusters, cell types, etc). If NULL (default), hdWGCNA uses the Seurat Idents as the group.
#' @param multi.group.by A string containing the name of a column in the Seurat object with groups for consensus WGCNA (dataset, sample, condition, etc)
#' @param multi_groups A character vecrtor containing the names of groups to select
#' @param assay The name of the assay in the Seurat object
#' @param slot The name of the slot in the Seurat object (counts, data)
#' @param mat A Matrix containing gene expression data. Supplying a matrix using this parameter ignores other options. This is almost exclusively used for pseudobulk analysis.
#' @param wgcna_name A string containing the name of the WGCNA slot in seurat_obj@misc. Default = NULL which retrieves the currently active WGCNA data
#' @keywords scRNA-seq
#' @export
SetMultiExpr <- function(
  seurat_obj,
  group_name,
  use_metacells=TRUE,
  group.by=NULL,
  multi.group.by = NULL,
  multi_groups = NULL,
  assay=NULL,
  slot='data',
  mat=NULL,
  mat_group_delim=3,
  wgcna_name=NULL,
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get the WGCNA genes:
  params <- GetWGCNAParams(seurat_obj, wgcna_name)
  gene_names <- GetWGCNAGenes(seurat_obj, wgcna_name)

  s_obj <- seurat_obj

  # get assay
  if(is.null(assay)){
    assay <- DefaultAssay(s_obj)
    warning(paste0('assay not specified, trying to use assay ', assay))
  }

  # get the different groups present if not specified by the user:
  if(is.null(multi_groups)){
    multi_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
  } else{
    seurat_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
    if(sum(multi_groups %in% seurat_groups) != length(multi_groups)){
      stop('Some or all groups specified in multi_groups not found in seurat_obj@meta.data[,multi.group.by]')
    }
  }

  # was a matrix supplied?
  if(is.null(mat)){

      # use metacells or whole seurat object?
    if(use_metacells){
      s_obj <- GetMetacellObject(seurat_obj, wgcna_name)
    } else{
      s_obj <- seurat_obj
    }

    # get the datExpr for each group
    datExpr_list <- lapply(multi_groups, function(cur_group){
      cur_datExpr <- SetDatExpr(
        seurat_obj,
        group_name = group_name,
        group.by = group.by,
        multi.group.by = multi.group.by,
        multi_group_name = cur_group,
        return_seurat = FALSE,
        use_metacells = use_metacells,
        wgcna_name = wgcna_name,
        assay = assay,
        slot = slot
      ) 
      as.matrix(cur_datExpr)
    })

  } else{
    sample_groups <- do.call(rbind, strsplit(rownames(mat), ':'))[,mat_group_delim]
    datExpr_list <- list()
    for(cur_group in multi_groups){
      cur_datExpr <- as.data.frame(mat[which(sample_groups == cur_group),])
      datExpr_list[[cur_group]]<- cur_datExpr
    }

  }

  # convert to multiExpr, get good genes:
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  genes_use <- WGCNA::goodGenesMS(multiExpr)
  gene_names <- gene_names[genes_use]

  # subset the multiExpr by the good genes::
  datExpr_list <- lapply(1:length(multiExpr), function(i){
    multiExpr[[i]]$data[,genes_use]
  })
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  names(multiExpr) <- multi_groups

  # update the WGCNA gene list:
  seurat_obj <- SetWGCNAGenes(seurat_obj, gene_names, wgcna_name)

  # set the multiExpr in the Seurat object
  seurat_obj@misc[[wgcna_name]]$multiExpr <- multiExpr
  seurat_obj

}


#' GetMultiExpr
#'
#' This function gets the expression matrix from the metacell object.
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
GetMultiExpr <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$multiExpr

}

############################
# WGCNA params
###########################

#' SetWGCNAParams
#'
#' @param seurat_obj A Seurat object
#' @param params list of WGCNA parameters
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetWGCNAParams <- function(seurat_obj, params, wgcna_name=NULL){

  if(is.null(GetActiveWGCNA(seurat_obj)$wgcna_params)){
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params <- params
  } else{
    for(i in 1:length(params)){
      seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_params[[names(params)[i]]] <- params[[i]]
    }
  }
  seurat_obj
}

#' GetWGCNAParams
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetWGCNAParams <- function(seurat_obj, wgcna_name=NULL){
  # test if wgcna_name is valid (TODO)

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_params
}

############################
# SoftPower Table
###########################

#' SetPowerTable
#'
#' @param seurat_obj A Seurat object
#' @param power_table a dataframe containing the results of the soft power test
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetPowerTable <- function(seurat_obj, power_table, wgcna_name=NULL){
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # add power table to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_powerTable <- power_table
  seurat_obj
}

#' GetPowerTable
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetPowerTable <- function(seurat_obj, wgcna_name=NULL){
  # get data from active assay if wgcna_name is not given

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_powerTable
}

############################
# WGCNA Network
###########################

#' SetNetworkData
#'
#' @param seurat_obj A Seurat object
#' @param net list of network data from WGCNA
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetNetworkData <- function(seurat_obj, net, wgcna_name=NULL){
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # add network data to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_net <- net
  seurat_obj
}


#' GetNetworkData
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetNetworkData <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_net
}

############################
# WGCNA modules dataframe
###########################


#' SetModules
#'
#' @param seurat_obj A Seurat object
#' @param modules dataframe containing gene module assignments
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetModules <- function(seurat_obj, mod_df, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # set module df
  seurat_obj@misc[[wgcna_name]]$wgcna_modules <- mod_df
  seurat_obj
}

#' GetModules
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetModules <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$wgcna_modules
}


#' GetHubGenes
#'
#' Extract the top N hub genes for a given set of modules. This function outputs
#' a table with the gene name, the module, and the kME for that module for the
#' top N hub genes.
#'
#' @param seurat_obj A Seurat object
#' @param n_hubs the number of hub genes to select for each module
#' @param mods list of modules, selects all modules by default
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetHubGenes <- function(
  seurat_obj,
  n_hubs = 10,
  mods = NULL,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 'grey')

  if(is.null(mods)){
    mods <- levels(modules$module); mods <- mods[mods != 'grey']
  } else{
    if(!all(mods %in% modules$module)){
      stop("Invalid selection for mods.")
    }
  }

  #get hub genes:
  hub_df <- do.call(rbind, lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur <- cur[,c('gene_name', 'module', paste0('kME_', cur_mod))]
    names(cur)[3] <- 'kME'
    cur <- dplyr::arrange(cur, kME)
    cur %>% dplyr::top_n(n_hubs, wt=kME)
  }))
  rownames(hub_df) <- 1:nrow(hub_df)
  hub_df

}


############################
# Module Eigengenes
###########################

#' SetMEs
#'
#' @param seurat_obj A Seurat object
#' @param MEs dataframe or matrix containing module eigengenes
#' @param harmonized logical indicating whether MEs have been harmonized
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMEs <- function(seurat_obj, MEs, harmonized=TRUE, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # harmonized MEs?
  if(harmonized){
    seurat_obj@misc[[wgcna_name]]$hMEs <- MEs
  } else{
    seurat_obj@misc[[wgcna_name]]$MEs <- MEs
  }
  seurat_obj
}

#' GetMEs
#'
#' Function to retrieve module eigengens from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param harmonized logical indicating whether MEs have been harmonized
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMEs <- function(seurat_obj, harmonized=TRUE, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get harmonized MEs?
  if(harmonized == TRUE && !is.null(seurat_obj@misc[[wgcna_name]]$hMEs)){
    MEs <- seurat_obj@misc[[wgcna_name]]$hMEs
  } else{
    MEs <- seurat_obj@misc[[wgcna_name]]$MEs
  }
  MEs
}


#' SetMELoadings
#'
#' @param seurat_obj A Seurat object
#' @param loadings named numeric vector with eigengene loadings
#' @param harmonized logical indicating whether MEs have been harmonized
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMELoadings <- function(seurat_obj, loadings, harmonized=TRUE, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # harmonized MEs?
  if(harmonized){
    seurat_obj@misc[[wgcna_name]]$hME_loadings <- c(seurat_obj@misc[[wgcna_name]]$hME_loadings, loadings)
  } else{
    seurat_obj@misc[[wgcna_name]]$ME_loadings <- c(seurat_obj@misc[[wgcna_name]]$ME_loadings, loadings)
  }
  seurat_obj
}

#' GetMELoadings
#'
#' Function to retrieve module eigengens from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param harmonized logical indicating whether MEs have been harmonized
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMELoadings <- function(seurat_obj, harmonized=TRUE, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get harmonized MEs?
  if(harmonized == TRUE && !is.null(seurat_obj@misc[[wgcna_name]]$hME_loadings)){
    MEs <- seurat_obj@misc[[wgcna_name]]$hME_loadings
  } else{
    MEs <- seurat_obj@misc[[wgcna_name]]$ME_loadings
  }
  MEs
}


############################
# GO term table
###########################

#' SetEnrichrTable
#'
#' @param seurat_obj A Seurat object
#' @param enrichr_table dataframe storing the results of running enrichr
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetEnrichrTable <- function(seurat_obj, enrich_table, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # set enrichr table
  seurat_obj@misc[[wgcna_name]]$enrichr_table <- enrich_table
  seurat_obj
}



#' GetEnrichrTable
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetEnrichrTable <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$enrichr_table
}


############################
# Module Scores
###########################


#' SetModuleScores
#'
#' @param seurat_obj A Seurat object
#' @param mod_scores dataframe storing the module expression scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetModuleScores <- function(seurat_obj, mod_scores, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$module_scores <- mod_scores
  seurat_obj
}

#' GetModuleScores
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetModuleScores <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$module_scores
}

############################
# Average Module Expression
###########################

#' SetAvgModuleExpr
#'
#' @param seurat_obj A Seurat object
#' @param avg_mods dataframe storing the average expression of all genes in the same module
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetAvgModuleExpr <- function(seurat_obj, avg_mods, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$avg_modules <- avg_mods
  seurat_obj
}

#' GetAvgModuleExpr
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetAvgModuleExpr <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$avg_modules
}

############################
# TOM
###########################

#' GetTOM
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTOM <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get WGCNA genes:
  gene_names <- GetWGCNAGenes(seurat_obj, wgcna_name)

  # load TOM
  tom_files <- GetNetworkData(seurat_obj, wgcna_name)$TOMFiles

  if(!file.exists(tom_files[[1]])){
    stop(paste0("TOM file ", tom_files[[1]], ' not found. Please update path to TOM file.'))
  }

  load(tom_files[[1]])

  TOM <- as.matrix(consTomDS)
  rownames(TOM) <- gene_names; colnames(TOM) <- gene_names
  TOM

}


############################
# TF Match Matrix (not stored within the WGCNA slot)
###########################

#' SetMotifMatrix
#'
#' @param seurat_obj A Seurat object
#' @param tf_match matrix containing tf-promoter matches
#' @keywords scRNA-seq
#' @export
#' @examples SetMotifMatrix
SetMotifMatrix <- function(seurat_obj, tf_match){

  # make a spot for the motif info if it's not already there:
  if(is.null(seurat_obj@misc$motifs)){
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$tf_match_matrix <- tf_match
  seurat_obj
}

#' GetMotifMatrix
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples GetMotifMatrix
GetMotifMatrix <- function(seurat_obj){
  seurat_obj@misc$motifs$tf_match_matrix
}


############################
# Motif table
###########################

#' SetMotifs
#'
#' @param seurat_obj A Seurat object
#' @param motif_df dataframe containing info about the motifs being analyzed
#' @keywords scRNA-seq
#' @export
#' @examples SetMotifs
SetMotifs <- function(seurat_obj, motif_df){

  # make a spot for the motif info if it's not already there:
  if(is.null(seurat_obj@misc$motifs)){
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$motif_df <- motif_df
  seurat_obj
}



#' GetMotifs
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples GetMotifs
GetMotifs <- function(seurat_obj){
  seurat_obj@misc$motifs$motif_df
}

############################
# PFM List
###########################


#' SetPFMList
#'
#' @param seurat_obj A Seurat object
#' @param pfm_list list of pfm objects
#' @keywords scRNA-seq
#' @export
#' @examples SetPFMList
SetPFMList <- function(seurat_obj, pfm_list){

  # make a spot for the motif info if it's not already there:
  if(is.null(seurat_obj@misc$motifs)){
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$pfm_list <- pfm_list
  seurat_obj
}

#' GetPFMList
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples GetPFMList
GetPFMList <- function(seurat_obj){
  seurat_obj@misc$motifs$pfm_list
}


############################
# TF Target Genes:
###########################

#' SetMotifTargets
#'
#' @param seurat_obj A Seurat object
#' @param motif_targets list of motifs and their target genes
#' @keywords scRNA-seq
#' @export
#' @examples SetMotifTargets
SetMotifTargets <- function(seurat_obj, motif_targets){

  # make a spot for the motif info if it's not already there:
  if(is.null(seurat_obj@misc$motifs)){
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$motif_targets <- motif_targets
  seurat_obj
}


#' GetMotifTargets
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples GetMotifTargets
GetMotifTargets <- function(seurat_obj){
  seurat_obj@misc$motifs$motif_targets
}


############################
# motif overlap
###########################


#' SetMotifOverlap
#'
#' @param seurat_obj A Seurat object
#' @param overlap_df dataframe containing motif-module overlap info
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMotifOverlap <- function(seurat_obj, overlap_df, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$motif_module_overlaps <- overlap_df
  seurat_obj
}


#' GetMotifOverlap
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMotifOverlap <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$motif_module_overlaps
}


############################
# motif scores
###########################


#' SetMotifScores
#'
#' @param seurat_obj A Seurat object
#' @param tf_scores dataframe of tf motif target scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMotifScores <- function(seurat_obj, tf_scores, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$motif_target_scores <- tf_scores
  seurat_obj
}


#' GetMotifScores
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMotifScores <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$motif_target_scores
}

############################
# ModuleUMAP
###########################

#' SetModuleUMAP
#'
#' @param seurat_obj A Seurat object
#' @param umap_df dataframe of UMAP coordinates
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetModuleUMAP <- function(seurat_obj, umap_df, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$module_umap <- umap_df
  seurat_obj
}

#' GetModuleUMAP
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetModuleUMAP <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$module_umap
}

############################
# ModuleTraitCorrelation
###########################


#' SetModuleTraitCorrelation
#'
#' @param seurat_obj A Seurat object
#' @param mt_cor matrix of module-trait correlation results
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetModuleTraitCorrelation <- function(seurat_obj, mt_cor, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$mt_cor <- mt_cor
  seurat_obj
}

#' GetModuleTraitCorrelation
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetModuleTraitCorrelation <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$mt_cor
}

############################
# ModulePreservation
###########################


#' SetModulePreservation
#'
#' @param seurat_obj A Seurat object
#' @param mt_cor matrix of module-trait correlation results
#' @param mod_name name of the module preservation test to store
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetModulePreservation <- function(seurat_obj, mod_pres, mod_name, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # make an empty list if module preservation hasn't been called yet
  if(is.null(seurat_obj@misc[[wgcna_name]]$module_preservation)){
    seurat_obj@misc[[wgcna_name]]$module_preservation <- list()
  }

  seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]] <- mod_pres
  seurat_obj
}



#' GetModulePreservation
#'
#' @param seurat_obj A Seurat object
#' @param mod_name name of the module preservation test to store
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetModulePreservation <- function(seurat_obj, mod_name, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  if(is.null(seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]])){
    stop("Invalid module preservation name.")
  }
  seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]]
}


############################
# Reset module names:
###########################


#' ResetModuleNames
#'
#' Reset the uname of each hdWGCNA module
#'
#' @param seurat_obj A Seurat object
#' @param new_name string containing the base name to re-name the modules
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
ResetModuleNames <- function(
  seurat_obj,
  new_name = "M",
  reset_levels = FALSE,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  old_mods <- levels(modules$module)
  if('grey' %in% modules$module){
    nmods <- length(old_mods) - 1
  } else{
      nmods <- length(old_mods)
  }

  # if it's a named list:
  if(class(new_name) == 'list'){
    if(all(names(new_name) %in% old_mods)){
      ix <- match(names(new_name), old_mods)
      new_names <- old_mods
      new_names[ix] <- as.character(new_name)
      new_names <- new_names[new_names != 'grey']
    } else{
      stop("Some module names present in new_name are not found in this hdWGCNA experiment.")
    }
  } else if(length(new_name) == 1){
    new_names <- paste0(new_name, 1:nmods)
  } else if(length(new_name) == nmods){
    new_names <- new_name
  } else{
    stop("Invalid input for new_name.")
  }

  if('grey' %in% modules$module){

    grey_ind <- which(old_mods == 'grey')

    # account for when grey is first / last
    if(grey_ind == 1){
      new_names <- c('grey', new_names)
    } else if(grey_ind == length(old_mods)){
      new_names <- c(new_names, 'grey')
    } else{
      new_names <- c(new_names[1:(grey_ind-1)], 'grey', new_names[grey_ind:length(new_names)])
    }

  }

  # update kMEs
  new_kMEs <- paste0('kME_', new_names)
  colnames(modules) <- c(colnames(modules)[1:3], new_kMEs)

  # update module names
  new_mod_df <- data.frame(
    old = old_mods ,
    new = new_names
  )

  modules$module <- factor(
    new_mod_df[match(modules$module, new_mod_df$old),'new'],
    levels = as.character(new_mod_df$new)
  )

  # set module table
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

  # update hME table:
  hMEs <- GetMEs(seurat_obj, harmonized=TRUE, wgcna_name)
  if(!is.null(hMEs)){
    me_colnames <- colnames(hMEs)
    ix <- match(me_colnames, new_mod_df$old)
    colnames(hMEs) <- new_mod_df$new[ix]
    seurat_obj <- SetMEs(seurat_obj, hMEs, harmonized=TRUE, wgcna_name)
  }

  # update ME table
  MEs <- GetMEs(seurat_obj, harmonized=FALSE, wgcna_name)
  if(!is.null(MEs)){
    me_colnames <- colnames(MEs)
    ix <- match(me_colnames, new_mod_df$old)
    colnames(MEs) <- new_mod_df$new[ix]
    seurat_obj <- SetMEs(seurat_obj, MEs, harmonized=FALSE, wgcna_name)
  }

  # update module scores:
  module_scores <- GetModuleScores(seurat_obj, wgcna_name)
  if(!is.null(module_scores)){
    if(!("grey" %in% colnames(module_scores))){
      colnames(module_scores) <- new_mod_df$new[new_mod_df$new != 'grey']
    } else {
      colnames(module_scores) <- new_mod_df$new
    }
    seurat_obj <- SetModuleScores(seurat_obj, module_scores, wgcna_name)
  }

  # update average module expression:
  avg_exp <- GetAvgModuleExpr(seurat_obj, wgcna_name)
  if(!is.null(avg_exp)){
    if(!("grey" %in% colnames(avg_exp))){
      colnames(avg_exp) <- new_mod_df$new[new_mod_df$new != 'grey']
    } else {
      colnames(avg_exp) <- new_mod_df$new
    }
    seurat_obj <- SetAvgModuleExpr(seurat_obj, avg_exp, wgcna_name)
  }

  # update enrichr table:
  enrich_table <- GetEnrichrTable(seurat_obj, wgcna_name)
  if(!is.null(enrich_table)){
    enrich_table$module <- factor(
      new_mod_df[match(enrich_table$module, new_mod_df$old),'new'],
      levels = as.character(new_mod_df$new)
    )
    seurat_obj <- SetEnrichrTable(seurat_obj, enrich_table, wgcna_name)
  }

  # update motif overlap
  overlap_df <- GetMotifOverlap(seurat_obj, wgcna_name)
  if(!is.null(overlap_df)){
    overlap_df$module <- factor(
      new_mod_df[match(overlap_df$module, new_mod_df$old),'new'],
      levels = as.character(new_mod_df$new)
    )
    seurat_obj <- SetMotifOverlap(seurat_obj, overlap_df, wgcna_name)
  }

  # update module umap:
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  if(!is.null(umap_df)){
    umap_df$module <- factor(
      new_mod_df[match(umap_df$module, new_mod_df$old),'new'],
      levels = as.character(new_mod_df$new)
    )
    seurat_obj <- SetModuleUMAP(seurat_obj, umap_df, wgcna_name)
  }


  seurat_obj

}


############################
# Reset module names:
###########################

#' ResetModuleColors
#'
#' Reset the unique color for each hdWGCNA module
#'
#' @param seurat_obj A Seurat object
#' @param new_colors a character vector containing the new colors
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
ResetModuleColors <- function(
  seurat_obj,
  new_colors,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  mod_colors_df <- dplyr::select(modules, c(module, color)) %>%
    distinct %>% arrange(module)
  mod_colors <- mod_colors_df$color
  if('grey' %in% modules$mod){
    grey_ind <- which(mod_colors == 'grey')
  } else{
    grey_ind <- NA
  }

  # case where we give a named list:
  if(class(new_colors) == 'list'){
    ix <- match(names(new_colors), mod_colors_df$module)
    mod_colors_df[ix, 'color'] <- as.character(new_colors)
    new_colors <- mod_colors_df$color
  } else{
    if(is.na(grey_ind)){
      new_colors <- new_colors
    } else if(grey_ind == 1){
      new_colors <- c('grey', new_colors)
    } else if(grey_ind == length(mod_colors)){
      new_colors <- c(new_colors, 'grey')
    } else{
      new_colors <- c(new_colors[1:(grey_ind-1)], 'grey', new_colors[grey_ind:length(new_colors)])
    }
  }

  new_color_df <- data.frame(
    old = mod_colors,
    new = new_colors
  )

  modules$color <- new_color_df[match(modules$color, new_color_df$old),'new']

  # set module table
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)

  # update motif overlap
  overlap_df <- GetMotifOverlap(seurat_obj, wgcna_name)
  if(!is.null(overlap_df)){
    overlap_df$color <- new_color_df[match(overlap_df$color, new_color_df$old),'new']
    seurat_obj <- SetMotifOverlap(seurat_obj, overlap_df, wgcna_name)
  }

  # update module umap:
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  if(!is.null(umap_df)){
    umap_df$color <- new_color_df[match(umap_df$color, new_color_df$old),'new']
    seurat_obj <- SetModuleUMAP(seurat_obj, umap_df, wgcna_name)
  }

  # update hdWGCNA dendrogram:
  net <- GetNetworkData(seurat_obj, wgcna_name)
  net$colors <- new_color_df[match(net$colors, new_color_df$old),'new']
  seurat_obj <- SetNetworkData(seurat_obj, net, wgcna_name)

  seurat_obj

}
