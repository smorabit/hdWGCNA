
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
  check <- wgcna_name %in% names(seurat_obj@misc) 
  if(!check){
    stop(paste0("Invalid wgcna_name supplied: ", wgcna_name))
  }  
}


# # get any WGCNA data, but by default get the active
# GetWGCNA <- function(seurat_obj, wgcna_name=NULL){

#   # test if wgcna_name is valid (TODO)

#   # get data from active assay if wgcna_name is not given
#   if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

#   seurat_obj@misc[[wgcna_name]]
# }

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)
  
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

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_genes
}


#' Set Expression Data (Standard WGCNA)
#'
#' This function sets up the expression matrix input for standard (non-consensus) WGCNA
#' based on the metacell expression matrix, the full expression matrix, or a provided
#' pseudobulk expression matrix.
#'
#' @description
#' The function operates in three modes:
#' 1. **Internal Seurat/Metacell Mode (Default):** Extracts expression data directly from the
#'    Seurat object (either single-cell or metacell data).
#' 2. **Pseudobulk Mode (SummarizedExperiment):** Extracts expression data from a provided
#'    SummarizedExperiment object (passed to `mat`). This is the recommended approach for
#'    pseudobulk analysis.
#' 3. **External Matrix Mode:** Sets the expression data from a provided matrix (passed to `mat`).
#'
#' @param seurat_obj A Seurat object containing the hdWGCNA experiment.
#' @param group_name A string containing the group to subset the data by (e.g., a specific
#'   cluster or cell type). Only used if pulling data from `seurat_obj`.
#' @param use_metacells Logical; if TRUE (default), use the metacell expression matrix.
#'   If FALSE, use the full single-cell expression matrix. Ignored if `mat` is provided.
#' @param group.by A string containing the name of a column in the Seurat object with
#'   cell groups (clusters, cell types, etc). If NULL (default), uses Seurat Idents.
#' @param multi.group.by A string containing the name of a column in the Seurat object
#'   with groups to subset further (e.g. dataset, sample). Only used if pulling data
#'   from `seurat_obj`.
#' @param multi_group_name A string or character vector specifying which groups from
#'   `multi.group.by` to include.
#' @param return_seurat Logical; if TRUE (default), returns the Seurat object with the
#'   `datExpr` slot populated. If FALSE, returns the data frame of expression data.
#' @param assay The name of the assay in the Seurat object (e.g., "RNA", "SCT").
#' @param slot The name of the slot in the Seurat object (e.g., "counts", "data").
#'   Used for Seurat v4 compatibility.
#' @param layer The name of the layer in the Seurat object (e.g., "counts", "data")
#'   OR the name of the assay in the SummarizedExperiment (e.g., "VST", "counts")
#'   if `mat` is provided.
#' @param mat A Matrix or SummarizedExperiment object containing gene expression data.
#'   - **SummarizedExperiment:** The function extracts the assay specified by `layer`
#'     and subsets genes to match `GetWGCNAGenes(seurat_obj)`.
#'   - **Matrix:** The function assumes columns are genes and rows are samples.
#' @param features A character vector of genes to use. If NULL (default), uses the
#'   genes stored in the hdWGCNA experiment.
#' @param wgcna_name A string containing the name of the WGCNA slot in `seurat_obj@misc`.
#'   Default = NULL, which retrieves the currently active WGCNA data.
#' @param ... Additional arguments passed to `WGCNA::goodGenes`.
#'
#' @return A Seurat object with the `datExpr` slot populated in the specified `wgcna_name`
#'   experiment, or a data frame if `return_seurat = FALSE`.
#'
#' @details
#' This function automatically runs `WGCNA::goodGenes` to exclude genes with zero variance
#' or excessive missingness. If `mat` is a SummarizedExperiment, the function automatically
#' transposes the assay to the required (Samples x Genes) format.
#'
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
  layer = 'data',
  mat=NULL,
  features=NULL,
  wgcna_name=NULL,
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # Check assay validity (only strictly needed if pulling from Seurat)
  if(is.null(assay)){
      assay <- DefaultAssay(seurat_obj)
      if(is.null(mat)) warning(paste0('assay not specified, trying to use assay ', assay))
  }

  # get parameters from seurat object
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  if(is.null(features)){
    genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  } else{
    if(all(features %in% rownames(seurat_obj))){
      genes_use <- features 
    } else{
      stop('Some features not found in rownames(seurat_obj).')
    }
  }

  # -------------------------------------------------------
  # Logic: Determine Source of Data
  # -------------------------------------------------------

  # Case 1: Pulling from Seurat/Metacells (No external matrix)
  if(is.null(mat)){
      
    # check that selected assay is in the seurat object 
    if(!(assay %in% names(seurat_obj))){
      stop(paste0('Invalid choice of assay: ', assay, ' not found in Assays(seurat_obj).'))
    }
    
    # check that slot is valid 
    if(!(slot %in% c('counts', 'data', 'scale.data'))){
      stop('Invalid choice of slot. Valid choices are counts, data, or scale.data.')
    }

    # check that layer is valid 
    if(!(layer %in% c('counts', 'data', 'scale.data'))){
      stop('Invalid choice of layer. For seurat objects, valid choices are counts, data, or scale.data.')
    }

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
      if(!(group.by %in% colnames(s_obj@meta.data))){
        m_cell_message <- ""
        if(use_metacells){m_cell_message <- "metacell"}
        stop(paste0(group.by, ' not found in the meta data of the ', m_cell_message, ' Seurat object'))
      }
      if(!all(group_name %in% s_obj@meta.data[[group.by]])){
        groups_not_found <- group_name[!(group_name %in% s_obj@meta.data[[group.by]])]
        stop(paste0("Some groups in group_name are not found in the seurat_obj: ", paste(groups_not_found, collapse=', ')))
      }
      seurat_meta <- seurat_meta %>% subset(get(group.by) %in% group_name)
    }

    # subset further if multiExpr:
    if(!is.null(multi.group.by)){
      if(!(multi.group.by %in% colnames(s_obj@meta.data))){
        stop(paste0(multi.group.by, ' not found in the meta data.'))
      }
      seurat_meta <- seurat_meta %>% subset(get(multi.group.by) %in% multi_group_name)
    }

    # get list of cells to use
    cells <- rownames(seurat_meta)

    # get expression data from seurat obj
    if(CheckSeurat5()){
      exp <- SeuratObject::LayerData(s_obj, assay=assay, layer=layer)
    } else{
      exp <- Seurat::GetAssayData(s_obj, assay=assay, slot=slot)
    }
    
    # Subset to WGCNA genes and selected cells
    # We transpose here to get Samples x Genes
    datExpr <- as.data.frame(t(as.matrix(exp[genes_use, cells, drop=FALSE])))
    
  } else {
      
    # Case 2: SummarizedExperiment provided
    if(inherits(mat, "SummarizedExperiment")){
        
        # Check if the requested assay/layer exists
        # We reuse the 'layer' argument here to select the SE assay (e.g. 'VST')
        if(!(layer %in% SummarizedExperiment::assayNames(mat))){
            stop(paste0("Assay '", layer, "' not found in SummarizedExperiment. Available assays: ", 
                        paste(SummarizedExperiment::assayNames(mat), collapse=", ")))
        }
        
        # Intersect SE genes with WGCNA genes
        genes_keep <- intersect(rownames(mat), genes_use)
        
        if(length(genes_keep) == 0){
             stop("No intersection between SummarizedExperiment rownames and selected WGCNA genes.")
        }
        
        # Extract, subset genes, and transpose to (Samples x Genes)
        # Note: We rely on the user to have filtered samples/pseudobulks in the SE object 
        # prior to calling this function if they wanted to subset groups.
        datExpr <- t(as.matrix(SummarizedExperiment::assay(mat, layer)[genes_keep, , drop=FALSE]))
        datExpr <- as.data.frame(datExpr)
        
    } else {
        
        # Case 3: Standard Matrix provided
        datExpr <- mat

        # cast it to a dataframe
        if(!is.data.frame(datExpr)){
          datExpr <- as.data.frame(datExpr)
        }

        # are the colnames genes?
        if(!all(colnames(datExpr) %in% rownames(seurat_obj))){
          stop("colnames of the provided matrix are invalid. Make sure that the colnames are features (genes), and that all of these features are in the seurat_obj")
        }
    }
  }

  if(return_seurat){
    
    # Run WGCNA's goodGenes check
    # We pass '...' so users can control verbose, minFraction, etc.
    is_good <- WGCNA::goodGenes(datExpr, ...)
    
    if(sum(is_good) < 2) {
        stop("Too few genes remaining after goodGenes check.")
    }
    
    # Subset datExpr
    datExpr <- datExpr[, is_good]
    
    # Update the gene list to match valid genes
    gene_list <- colnames(datExpr)

    # update the WGCNA gene list in the object:
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
  CheckWGCNAName(seurat_obj, wgcna_name)
  
  seurat_obj@misc[[wgcna_name]]$datExpr

}


' Set Multi-Set Expression Data (Consensus WGCNA)
#'
#' This function prepares the expression data for Consensus WGCNA analysis. It populates
#' the `multiExpr` slot in the Seurat object, which contains a list of expression matrices
#' (one for each consensus group, e.g., dataset, sample, condition).
#'
#' @description
#' The function operates in three modes:
#' 1. **Internal Seurat/Metacell Mode (Default):** Extracts expression data directly from the
#'    Seurat object (either single-cell or metacell data) based on `group_name`.
#' 2. **Pseudobulk Mode (SummarizedExperiment):** Extracts expression data from a provided
#'    SummarizedExperiment object (passed to `mat`). This is the recommended approach for
#'    pseudobulk consensus analysis.
#' 3. **External Matrix Mode:** Extracts expression data from a provided large matrix
#'    (passed to `mat`) where row names contain delimited group identifiers.
#'
#' @param seurat_obj A Seurat object containing the hdWGCNA experiment.
#' @param group_name A string containing the specific group to analyze (e.g., a specific
#'   cluster or cell type). This filters the data when using Internal Seurat/Metacell Mode.
#' @param use_metacells Logical; if TRUE (default), use the metacell expression matrix
#'   stored in the hdWGCNA experiment. If FALSE, use the full single-cell expression matrix.
#'   Ignored if `mat` is provided.
#' @param group.by A string containing the name of a column in the Seurat object with
#'   cell groups (clusters, cell types, etc). If NULL (default), uses Seurat Idents.
#' @param multi.group.by A string containing the name of the column that defines the
#'   consensus groups (e.g., "dataset", "sample", "condition").
#'   - If using **Seurat/Metacells**, this must be a column in `seurat_obj@meta.data`.
#'   - If using **Pseudobulk (SE)**, this must be a column in `colData(mat)`.
#' @param multi_groups A character vector specifying which groups from `multi.group.by`
#'   to include. If NULL, all unique groups are used.
#' @param assay The name of the assay in the Seurat object (e.g., "RNA", "SCT").
#' @param slot The name of the slot in the Seurat object (e.g., "counts", "data").
#'   Used for Seurat v4 compatibility.
#' @param layer The name of the layer in the Seurat object (e.g., "counts", "data")
#'   OR the name of the assay in the SummarizedExperiment (e.g., "VST", "counts")
#'   if `mat` is provided.
#' @param mat A Matrix or SummarizedExperiment object containing gene expression data.
#'   - **SummarizedExperiment:** The function extracts the assay specified by `layer`
#'     and subsets columns based on `multi.group.by` in `colData`.
#'   - **Matrix:** The function assumes row names are delimited (e.g., "Cluster1:SampleA")
#'     and splits them using `mat_group_delim`.
#' @param mat_group_delim Character; the delimiter used in the row names of `mat`
#'   if `mat` is a matrix (default is ":"). Ignored if `mat` is a SummarizedExperiment.
#' @param wgcna_name A string containing the name of the WGCNA slot in `seurat_obj@misc`.
#'   Default = NULL, which retrieves the currently active WGCNA data.
#' @param ... Additional arguments passed to internal helper functions.
#'
#' @return A Seurat object with the `multiExpr` slot populated in the specified `wgcna_name`
#'   experiment.
#'
#' @details
#' This function automatically aligns the genes in the provided data (`mat` or `seurat_obj`)
#' with the genes selected for WGCNA (via `SetupForWGCNA` or `GetWGCNAGenes`). It also
#' runs `WGCNA::goodGenesMS` to exclude genes with zero variance or excessive missingness
#' across the consensus groups.
#'
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
  layer = 'data',
  mat=NULL,
  mat_group_delim=3,
  wgcna_name=NULL,
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

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
    if(is.null(mat) || !inherits(mat, "SummarizedExperiment")){
       multi_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
    }
  } else{
    # Validate groups if we are using the internal Seurat object
    if(is.null(mat)){
        seurat_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
        if(sum(multi_groups %in% seurat_groups) != length(multi_groups)){
          stop('Some or all groups specified in multi_groups not found in seurat_obj@meta.data[,multi.group.by]')
        }
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
        slot = slot,
        layer = layer
      ) 
      as.matrix(cur_datExpr)
    })

  } else {
      
    # CASE 2: SummarizedExperiment provided (Pseudobulk)
    if(inherits(mat, "SummarizedExperiment")){
        
        # Validation
        if(is.null(multi.group.by)){
            stop("You must provide 'multi.group.by' (the column name in colData) when using a SummarizedExperiment.")
        }
        if(!(multi.group.by %in% names(SummarizedExperiment::colData(mat)))){
            stop(paste0("Column '", multi.group.by, "' not found in SummarizedExperiment colData."))
        }
        
        # Check if the requested assay/layer exists
        # NOTE: We use the 'layer' argument to select the SE assay (e.g., 'VST', 'counts')
        if(!(layer %in% SummarizedExperiment::assayNames(mat))){
            stop(paste0("Assay '", layer, "' not found in SummarizedExperiment. Available assays: ", 
                        paste(SummarizedExperiment::assayNames(mat), collapse=", ")))
        }

        # Determine groups if not provided
        group_vec <- SummarizedExperiment::colData(mat)[[multi.group.by]]
        if(is.null(multi_groups)){
            multi_groups <- unique(as.character(group_vec))
        }

        # Extract data list
        datExpr_list <- lapply(multi_groups, function(x){
            
            # Identify columns for this group
            cells_keep <- group_vec == x
            
            # Subset the SE object to the WGCNA genes and the group cells
            # We use 'gene_names' here to ensure alignment with SetupForWGCNA
            genes_keep <- intersect(rownames(mat), gene_names)
            
            if(length(genes_keep) == 0){
                stop("No intersection between SummarizedExperiment rownames and selected WGCNA genes.")
            }
            
            # Extract, subset, and transpose to (Samples x Genes)
            dat <- SummarizedExperiment::assay(mat, layer)[genes_keep, cells_keep, drop=FALSE]
            t(as.matrix(dat))
        })
        names(datExpr_list) <- multi_groups

    # CASE 3: Large Matrix provided (Delimited rownames)
    } else {
        
        # Ensure we are splitting by the correct delimiter
        sample_groups <- do.call(rbind, strsplit(rownames(mat), ':'))[,mat_group_delim]
        datExpr_list <- list()
        
        for(cur_group in multi_groups){
          cur_datExpr <- as.data.frame(mat[which(sample_groups == cur_group),])
          datExpr_list[[cur_group]]<- cur_datExpr
        }
    }

  }

  # convert to multiExpr, get good genes:
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  genes_use <- WGCNA::goodGenesMS(multiExpr)
  
  # Update the gene_names based on the goodGenes check
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
  
  return(seurat_obj)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$multiExpr

}

############################
# WGCNA params
###########################


#' SetMetacellParams
#'
#' @param seurat_obj A Seurat object
#' @param params list of WGCNA parameters
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMetacellParams <- function(seurat_obj, params, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$metacell_params <- params
  seurat_obj
}

#' GetMetacellParams
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMetacellParams <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$metacell_params
}

#' SetWGCNAParams
#'
#' @param seurat_obj A Seurat object
#' @param params list of WGCNA parameters
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetWGCNAParams <- function(seurat_obj, params, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_params <- params
  seurat_obj
}

#' GetWGCNAParams
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetWGCNAParams <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  
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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)
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
SetModules <- function(seurat_obj, modules, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set module df
  seurat_obj@misc[[wgcna_name]]$wgcna_modules <- modules
  seurat_obj
}

#' GetModules
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
GetModules <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_modules
}





#' SetDegrees
#'
#' @param seurat_obj A Seurat object
#' @param degree_df dataframe containing gene module assignments
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
SetDegrees <- function(seurat_obj, degree_df, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set module df
  seurat_obj@misc[[wgcna_name]]$wgcna_degrees <- degree_df
  seurat_obj
}

#' GetDegrees
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
GetDegrees <- function(seurat_obj, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$wgcna_degrees
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
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get the modules table
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
    cur <- dplyr::arrange(cur, desc(kME))
    cur %>% dplyr::slice_max(n=n_hubs, order_by=kME)
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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
#' @param enrich_table dataframe storing the results of running enrichr
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetEnrichrTable <- function(seurat_obj, enrich_table, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$enrichr_table
}


#' SetEnrichRegulonTable
#'
#' @param seurat_obj A Seurat object
#' @param enrich_table dataframe storing the results of running enrichr
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetEnrichrRegulonTable <- function(seurat_obj, enrich_table, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set enrichr table
  seurat_obj@misc[[wgcna_name]]$enrichr_regulon_table <- enrich_table
  seurat_obj
}



#' GetEnrichrRegulonTable
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetEnrichrRegulonTable <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$enrichr_regulon_table
}



############################
# Module Scores
###########################


#' SetModuleScores
#'
#' @param seurat_obj A Seurat object
#' @param mod_scores dataframe storing the module expression scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
SetModuleScores <- function(seurat_obj, mod_scores, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$module_scores <- mod_scores
  seurat_obj
}

#' GetModuleScores
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
GetModuleScores <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

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
#' @export
SetAvgModuleExpr <- function(seurat_obj, avg_mods, wgcna_name=NULL){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get modules 
  modules <- GetModules(seurat_obj, wgcna_name)
  gene_names <- modules$gene_name

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

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
  CheckWGCNAName(seurat_obj, wgcna_name)

  if(is.null(seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]])){
    stop("Invalid module preservation name.")
  }
  seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]]
}




#' SetRegulonScores
#'
#' @param seurat_obj A Seurat object
#' @param regulon_scores dataframe storing the TF regulon scores 
#' @param target_type dataframe storing the TF regulon scores 
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetRegulonScores <- function(seurat_obj, regulon_scores, target_type, wgcna_name=NULL){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)


    # if regulon scores have not been set, make a list to store them
    if(is.null(seurat_obj@misc[[wgcna_name]]$regulon_scores)){
        tmp <- list(regulon_scores); names(tmp) <- target_type
        seurat_obj@misc[[wgcna_name]]$regulon_scores <- tmp
    } else{
        seurat_obj@misc[[wgcna_name]]$regulon_scores[[target_type]] <- regulon_scores
    }
    seurat_obj
}

#' GetRegulonScores
#'
#' @param seurat_obj A Seurat object
#' @param target_type dataframe storing the TF regulon scores 
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetRegulonScores <- function(seurat_obj, target_type, wgcna_name=NULL){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)
    
    # get the regulon scores
    seurat_obj@misc[[wgcna_name]]$regulon_scores[[target_type]] 
}


############################
# getters and setters
###########################

#' SetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param tf_net dataframe storing the TF network info in ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFNetwork <- function(seurat_obj, tf_net, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_net <- tf_net
  seurat_obj
}

#' GetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFNetwork <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_net 
}


#' SetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param tf_eval dataframe storing the TF network evaluation info from ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFEval <- function(seurat_obj, tf_eval, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_eval <- tf_eval
  seurat_obj
}

#' GetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFEval <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_eval
}

#' SetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param tf_regulons dataframe storing the TF regulon info from AssignTFRegulons
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFRegulons <- function(seurat_obj, tf_regulons, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_regulons <- tf_regulons
  seurat_obj
}

#' GetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFRegulons <- function(seurat_obj, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_regulons

}


############################
# Reset module names:
###########################
