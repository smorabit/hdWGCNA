############################
# Active WGCNA
###########################

#' SetActiveWGCNA
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment to set as active. An empty experiment is created in seurat_obj@misc if it does not exist yet.
#' @return The input Seurat object with \code{seurat_obj@misc$active_wgcna} set to \code{wgcna_name}.
#' @keywords scRNA-seq
#' @export
SetActiveWGCNA <- function(seurat_obj, wgcna_name) {
  # set the active_wgcna variable
  seurat_obj@misc$active_wgcna <- wgcna_name

  # initialize empty list for this WGCNA if it doesn't exist yet
  if (!(wgcna_name %in% names(seurat_obj@misc))) {
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]] <- list()
  }
  seurat_obj
}


#' GetActiveWGCNA
#'
#' @param seurat_obj A Seurat object.
#' @return A list containing all data stored for the active hdWGCNA experiment.
#' @keywords scRNA-seq
#' @export
GetActiveWGCNA <- function(seurat_obj) {
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]
}


#' GetActiveWGCNAName
#'
#' @param seurat_obj A Seurat object.
#' @return A character string with the name of the active hdWGCNA experiment, or NULL if none is set.
#' @keywords scRNA-seq
#' @export
GetActiveWGCNAName <- function(seurat_obj) {
  seurat_obj@misc$active_wgcna
}


#' CheckWGCNAName
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment to check for in the seurat_obj@misc slot.
#' @return Invisibly returns TRUE if \code{wgcna_name} is a valid hdWGCNA experiment; otherwise throws an error.
#' @keywords scRNA-seq
#' @keywords internal
#' @noRd
CheckWGCNAName <- function(seurat_obj, wgcna_name) {
  if (is.null(wgcna_name) || length(wgcna_name) != 1) {
    stop('No hdWGCNA experiment found. Run SetupForWGCNA() first, or provide wgcna_name.')
  }
  check <- wgcna_name %in% names(seurat_obj@misc)
  if (!check) {
    stop(paste0('Invalid wgcna_name supplied: ', wgcna_name))
  }
  invisible(TRUE)
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
#' @param seurat_obj A Seurat object.
#' @param metacell_obj A metacell Seurat object, or the name of another hdWGCNA experiment whose metacell object should be used.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the metacell object stored in \code{seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj}.
#' @keywords scRNA-seq
#' @export
SetMetacellObject <- function(seurat_obj, metacell_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # add metacell obj to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj <- metacell_obj
  seurat_obj
}


#' GetMetacellObject
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The metacell Seurat object for this experiment (following the reference if another experiment's name was stored), or NULL if no metacell object is found.
#' @keywords scRNA-seq
#' @export
GetMetacellObject <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  input_class <- class(seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj)
  if (input_class == 'Seurat') {
    return(seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj)
  } else if (input_class == 'character') {
    metacell_location <- seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj
    return(seurat_obj@misc[[metacell_location]]$wgcna_metacell_obj)
  } else {
    return(NULL)
  }
}


############################
# WGCNA genes
###########################

#' SetWGCNAGenes
#'
#' @param seurat_obj A Seurat object.
#' @param gene_list A character vector of genes to be used for WGCNA.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the gene list stored in \code{seurat_obj@misc[[wgcna_name]]$wgcna_genes}.
#' @keywords scRNA-seq
#' @export
SetWGCNAGenes <- function(seurat_obj, gene_list, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # add gene list to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_genes <- gene_list
  seurat_obj
}


#' GetWGCNAGenes
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A character vector of the genes selected for WGCNA.
#' @keywords scRNA-seq
#' @export
GetWGCNAGenes <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
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
  use_metacells = TRUE,
  group.by = NULL,
  multi.group.by = NULL,
  multi_group_name = NULL,
  return_seurat = TRUE,
  assay = NULL,
  slot = 'data',
  layer = 'data',
  mat = NULL,
  features = NULL,
  wgcna_name = NULL,
  ...
) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # Check assay validity (only strictly needed if pulling from Seurat)
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
    if (is.null(mat)) warning(paste0('assay not specified, trying to use assay ', assay))
  }

  # get parameters from seurat object
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  if (is.null(features)) {
    genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  } else {
    if (all(features %in% rownames(seurat_obj))) {
      genes_use <- features
    } else {
      stop('Some features not found in rownames(seurat_obj).')
    }
  }

  # -------------------------------------------------------
  # Logic: Determine Source of Data
  # -------------------------------------------------------

  # Case 1: Pulling from Seurat/Metacells (No external matrix)
  if (is.null(mat)) {
    # check that selected assay is in the seurat object
    if (!(assay %in% names(seurat_obj))) {
      stop(paste0('Invalid choice of assay: ', assay, ' not found in Assays(seurat_obj).'))
    }

    # check that slot is valid
    if (!(slot %in% c('counts', 'data', 'scale.data'))) {
      stop('Invalid choice of slot. Valid choices are counts, data, or scale.data.')
    }

    # check that layer is valid
    if (!(layer %in% c('counts', 'data', 'scale.data'))) {
      stop('Invalid choice of layer. For seurat objects, valid choices are counts, data, or scale.data.')
    }

    # get metacell object
    m_obj <- GetMetacellObject(seurat_obj, wgcna_name)

    # use metacells or whole seurat object?
    if (use_metacells && !is.null(m_obj)) {
      s_obj <- m_obj
    } else {
      if (is.null(m_obj)) {
        warning('Metacell Seurat object not found. Using full Seurat object instead.')
      }
      s_obj <- seurat_obj
    }

    # get the metadata from the seurat object:
    seurat_meta <- s_obj@meta.data

    # check the group.by params
    if (!is.null(group.by)) {
      if (!(group.by %in% colnames(s_obj@meta.data))) {
        m_cell_message <- ''
        if (use_metacells) {
          m_cell_message <- 'metacell'
        }
        stop(paste0(group.by, ' not found in the meta data of the ', m_cell_message, ' Seurat object'))
      }
      if (!all(group_name %in% s_obj@meta.data[[group.by]])) {
        groups_not_found <- group_name[!(group_name %in% s_obj@meta.data[[group.by]])]
        stop(paste0('Some groups in group_name are not found in the seurat_obj: ', paste(groups_not_found, collapse = ', ')))
      }
      seurat_meta <- seurat_meta %>% subset(get(group.by) %in% group_name)
    }

    # subset further if multiExpr:
    if (!is.null(multi.group.by)) {
      if (!(multi.group.by %in% colnames(s_obj@meta.data))) {
        stop(paste0(multi.group.by, ' not found in the meta data.'))
      }
      seurat_meta <- seurat_meta %>% subset(get(multi.group.by) %in% multi_group_name)
    }

    # get list of cells to use
    cells <- rownames(seurat_meta)

    # get expression data from seurat obj
    if (CheckSeurat5()) {
      exp <- SeuratObject::LayerData(s_obj, assay = assay, layer = layer)
    } else {
      exp <- Seurat::GetAssayData(s_obj, assay = assay, slot = slot)
    }

    # Subset to WGCNA genes and selected cells
    # We transpose here to get Samples x Genes
    datExpr <- as.data.frame(t(as.matrix(exp[genes_use, cells, drop = FALSE])))
  } else {
    # Case 2: SummarizedExperiment provided
    if (inherits(mat, 'SummarizedExperiment')) {
      # Check if the requested assay/layer exists
      # We reuse the 'layer' argument here to select the SE assay (e.g. 'VST')
      if (!(layer %in% SummarizedExperiment::assayNames(mat))) {
        stop(paste0(
          "Assay '", layer, "' not found in SummarizedExperiment. Available assays: ",
          paste(SummarizedExperiment::assayNames(mat), collapse = ', ')
        ))
      }

      # Intersect SE genes with WGCNA genes
      genes_keep <- intersect(rownames(mat), genes_use)

      if (length(genes_keep) == 0) {
        stop('No intersection between SummarizedExperiment rownames and selected WGCNA genes.')
      }

      # Extract, subset genes, and transpose to (Samples x Genes)
      # Note: We rely on the user to have filtered samples/pseudobulks in the SE object
      # prior to calling this function if they wanted to subset groups.
      datExpr <- t(as.matrix(SummarizedExperiment::assay(mat, layer)[genes_keep, , drop = FALSE]))
      datExpr <- as.data.frame(datExpr)
    } else {
      # Case 3: Standard Matrix provided
      datExpr <- mat

      # cast it to a dataframe
      if (!is.data.frame(datExpr)) {
        datExpr <- as.data.frame(datExpr)
      }

      # are the colnames genes?
      if (!all(colnames(datExpr) %in% rownames(seurat_obj))) {
        stop('colnames of the provided matrix are invalid. Make sure that the colnames are features (genes), and that all of these features are in the seurat_obj')
      }
    }
  }

  if (return_seurat) {
    # Run WGCNA's goodGenes check
    # We pass '...' so users can control verbose, minFraction, etc.
    is_good <- WGCNA::goodGenes(datExpr, ...)

    if (sum(is_good) < 2) {
      stop('Too few genes remaining after goodGenes check.')
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
  } else {
    out <- datExpr
  }
  out
}


#' GetDatExpr
#'
#' This function gets the WGCNA expression matrix.
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of the expression data used for WGCNA (samples or metacells by genes), as stored by \code{SetDatExpr}.
#' @keywords scRNA-seq
#' @export
GetDatExpr <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$datExpr
}


#' Set Multi-Set Expression Data (Consensus WGCNA)
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
  use_metacells = TRUE,
  group.by = NULL,
  multi.group.by = NULL,
  multi_groups = NULL,
  assay = NULL,
  slot = 'data',
  layer = 'data',
  mat = NULL,
  mat_group_delim = 3,
  wgcna_name = NULL,
  ...
) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get the WGCNA genes:
  params <- GetWGCNAParams(seurat_obj, wgcna_name)
  gene_names <- GetWGCNAGenes(seurat_obj, wgcna_name)

  s_obj <- seurat_obj

  # get assay
  if (is.null(assay)) {
    assay <- DefaultAssay(s_obj)
    warning(paste0('assay not specified, trying to use assay ', assay))
  }

  # get the different groups present if not specified by the user:
  if (is.null(multi_groups)) {
    if (is.null(mat) || !inherits(mat, 'SummarizedExperiment')) {
      multi_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
    }
  } else {
    # Validate groups if we are using the internal Seurat object
    if (is.null(mat)) {
      seurat_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
      if (sum(multi_groups %in% seurat_groups) != length(multi_groups)) {
        stop('Some or all groups specified in multi_groups not found in seurat_obj@meta.data[,multi.group.by]')
      }
    }
  }

  # was a matrix supplied?
  if (is.null(mat)) {
    # use metacells or whole seurat object?
    if (use_metacells) {
      s_obj <- GetMetacellObject(seurat_obj, wgcna_name)
    } else {
      s_obj <- seurat_obj
    }

    # get the datExpr for each group
    datExpr_list <- lapply(multi_groups, function(cur_group) {
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
    if (inherits(mat, 'SummarizedExperiment')) {
      # Validation
      if (is.null(multi.group.by)) {
        stop("You must provide 'multi.group.by' (the column name in colData) when using a SummarizedExperiment.")
      }
      if (!(multi.group.by %in% names(SummarizedExperiment::colData(mat)))) {
        stop(paste0("Column '", multi.group.by, "' not found in SummarizedExperiment colData."))
      }

      # Check if the requested assay/layer exists
      # NOTE: We use the 'layer' argument to select the SE assay (e.g., 'VST', 'counts')
      if (!(layer %in% SummarizedExperiment::assayNames(mat))) {
        stop(paste0(
          "Assay '", layer, "' not found in SummarizedExperiment. Available assays: ",
          paste(SummarizedExperiment::assayNames(mat), collapse = ', ')
        ))
      }

      # Determine groups if not provided
      group_vec <- SummarizedExperiment::colData(mat)[[multi.group.by]]
      if (is.null(multi_groups)) {
        multi_groups <- unique(as.character(group_vec))
      }

      # Extract data list
      datExpr_list <- lapply(multi_groups, function(x) {
        # Identify columns for this group
        cells_keep <- group_vec == x

        # Subset the SE object to the WGCNA genes and the group cells
        # We use 'gene_names' here to ensure alignment with SetupForWGCNA
        genes_keep <- intersect(rownames(mat), gene_names)

        if (length(genes_keep) == 0) {
          stop('No intersection between SummarizedExperiment rownames and selected WGCNA genes.')
        }

        # Extract, subset, and transpose to (Samples x Genes)
        dat <- SummarizedExperiment::assay(mat, layer)[genes_keep, cells_keep, drop = FALSE]
        t(as.matrix(dat))
      })
      names(datExpr_list) <- multi_groups

      # CASE 3: Large Matrix provided (Delimited rownames)
    } else {
      # Ensure we are splitting by the correct delimiter
      sample_groups <- do.call(rbind, strsplit(rownames(mat), ':'))[, mat_group_delim]
      datExpr_list <- list()

      for (cur_group in multi_groups) {
        cur_datExpr <- as.data.frame(mat[which(sample_groups == cur_group), ])
        datExpr_list[[cur_group]] <- cur_datExpr
      }
    }
  }

  # convert to multiExpr, get good genes:
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  genes_use <- WGCNA::goodGenesMS(multiExpr)

  # Update the gene_names based on the goodGenes check
  gene_names <- gene_names[genes_use]

  # subset the multiExpr by the good genes::
  datExpr_list <- lapply(seq_along(multiExpr), function(i) {
    multiExpr[[i]]$data[, genes_use]
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
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The multi-set expression data used for consensus WGCNA, as stored by \code{SetMultiExpr} (a WGCNA multiData list with one expression matrix per group).
#' @keywords scRNA-seq
#' @export
GetMultiExpr <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$multiExpr
}


############################
# WGCNA params
###########################


#' SetMetacellParams
#'
#' @param seurat_obj A Seurat object.
#' @param params A list of metacell construction parameters.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the metacell parameters stored in \code{seurat_obj@misc[[wgcna_name]]$metacell_params}.
#' @keywords scRNA-seq
#' @export
SetMetacellParams <- function(seurat_obj, params, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$metacell_params <- params
  seurat_obj
}


#' GetMetacellParams
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A list of the parameters used to construct metacells.
#' @keywords scRNA-seq
#' @export
GetMetacellParams <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$metacell_params
}


#' SetWGCNAParams
#'
#' @param seurat_obj A Seurat object.
#' @param params A list of WGCNA parameters.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the parameter list stored in \code{seurat_obj@misc[[wgcna_name]]$wgcna_params}.
#' @keywords scRNA-seq
#' @export
SetWGCNAParams <- function(seurat_obj, params, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_params <- params
  seurat_obj
}


#' GetWGCNAParams
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A list of the parameters stored for this hdWGCNA experiment.
#' @keywords scRNA-seq
#' @export
GetWGCNAParams <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_params
}


############################
# SoftPower Table
###########################

#' SetPowerTable
#'
#' @param seurat_obj A Seurat object.
#' @param power_table A data.frame containing the results of the soft power threshold test.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the power table stored in \code{seurat_obj@misc[[wgcna_name]]$wgcna_powerTable}.
#' @keywords scRNA-seq
#' @export
SetPowerTable <- function(seurat_obj, power_table, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # add power table to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_powerTable <- power_table
  seurat_obj
}


#' GetPowerTable
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of soft power threshold results, as computed by \code{TestSoftPowers}.
#' @keywords scRNA-seq
#' @export
GetPowerTable <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_powerTable
}


############################
# WGCNA Network
###########################

#' SetNetworkData
#'
#' @param seurat_obj A Seurat object.
#' @param net A list of network data returned by WGCNA network construction.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the network data stored in \code{seurat_obj@misc[[wgcna_name]]$wgcna_net}.
#' @keywords scRNA-seq
#' @export
SetNetworkData <- function(seurat_obj, net, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # add network data to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_net <- net
  seurat_obj
}


#' GetNetworkData
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A list of network data from \code{ConstructNetwork} (including module colors, dendrograms, and TOM file paths).
#' @keywords scRNA-seq
#' @export
GetNetworkData <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$wgcna_net
}


############################
# TOM
###########################

#' SetTOM
#'
#' @param seurat_obj A Seurat object.
#' @param TOM A gene-by-gene topological overlap matrix with gene names as row and column names.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the TOM stored in \code{seurat_obj@misc[[wgcna_name]]$wgcna_TOM}.
#' @keywords scRNA-seq
#' @export
SetTOM <- function(seurat_obj, TOM, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$wgcna_TOM <- TOM
  seurat_obj
}


#' GetTOM
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A gene-by-gene topological overlap matrix with gene names as row and column names. Uses the in-object TOM if present, otherwise loads it from the TOM file on disk.
#' @keywords scRNA-seq
#' @export
GetTOM <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # check if TOM is stored in-object first
  tom_in_obj <- seurat_obj@misc[[wgcna_name]]$wgcna_TOM
  if (!is.null(tom_in_obj)) {
    return(tom_in_obj)
  }

  # fall back to file-based loading
  modules <- GetModules(seurat_obj, wgcna_name)
  gene_names <- modules$gene_name

  tom_files <- GetNetworkData(seurat_obj, wgcna_name)$TOMFiles

  # objects made with hdWGCNA <= 0.4.13 and an absolute tom_outdir stored the path as
  # paste0(getwd(), '/', tom_outdir, ...), e.g. '/home/me//data/TOM/x_TOM.rda' or, on Windows,
  # 'C:/Users/me/D:/data/TOM/x_TOM.rda'; recover the absolute path at the end
  if (!file.exists(tom_files[[1]])) {
    fixed_path <- repair_tom_path(tom_files[[1]])
    if (!is.na(fixed_path) && file.exists(fixed_path)) {
      tom_files[[1]] <- fixed_path
    }
  }

  if (!file.exists(tom_files[[1]])) {
    stop(paste0('TOM file ', tom_files[[1]], ' not found. Please update path to TOM file.'))
  }

  load(tom_files[[1]])

  TOM <- as.matrix(consTomDS)
  rownames(TOM) <- gene_names
  colnames(TOM) <- gene_names
  TOM
}


############################
# WGCNA modules dataframe
###########################


#' SetModules
#'
#' @param seurat_obj A Seurat object.
#' @param modules A data.frame containing gene module assignments.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the module table stored in \code{seurat_obj@misc[[wgcna_name]]$wgcna_modules}.
#' @keywords scRNA-seq
#' @export
SetModules <- function(seurat_obj, modules, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set module df
  seurat_obj@misc[[wgcna_name]]$wgcna_modules <- modules
  seurat_obj
}


#' GetModules
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of gene module assignments, with columns for gene name, module, color, and kME values if computed.
#' @export
GetModules <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_modules
}


#' SetDegrees
#'
#' @param seurat_obj A Seurat object.
#' @param degree_df A data.frame containing intramodular connectivity (degree) for each gene.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the degree table stored in \code{seurat_obj@misc[[wgcna_name]]$wgcna_degrees}.
#' @export
SetDegrees <- function(seurat_obj, degree_df, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set module df
  seurat_obj@misc[[wgcna_name]]$wgcna_degrees <- degree_df
  seurat_obj
}


#' GetDegrees
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of intramodular connectivity (degree) for each gene, or NULL if not computed.
#' @export
GetDegrees <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$wgcna_degrees
}


#' GetHubGenes
#'
#' Extract the top N hub genes for a given set of modules. This function outputs
#' a table with the gene name, the module, and the kME for that module for the
#' top N hub genes.
#'
#' @param seurat_obj A Seurat object.
#' @param n_hubs The number of hub genes to select for each module.
#' @param mods A character vector of modules to include. Defaults to all modules except grey.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame with columns gene_name, module, and kME containing the top \code{n_hubs} genes by kME for each module.
#' @keywords scRNA-seq
#' @export
GetHubGenes <- function(
  seurat_obj,
  n_hubs = 10,
  mods = NULL,
  wgcna_name = NULL
) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get the modules table
  modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 'grey')

  if (is.null(mods)) {
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']
  } else {
    if (!all(mods %in% modules$module)) {
      stop('Invalid selection for mods.')
    }
  }

  # get hub genes:
  hub_df <- do.call(rbind, lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur <- cur[, c('gene_name', 'module', paste0('kME_', cur_mod))]
    names(cur)[3] <- 'kME'
    cur <- dplyr::arrange(cur, desc(kME))
    cur %>% dplyr::slice_max(n = n_hubs, order_by = kME)
  }))
  rownames(hub_df) <- seq_len(nrow(hub_df))
  hub_df
}


############################
# Module Eigengenes
###########################

#' SetMEs
#'
#' @param seurat_obj A Seurat object.
#' @param MEs A data.frame or matrix of module eigengenes (cells by modules).
#' @param harmonized Logical; if TRUE (default), store as harmonized MEs (hMEs), otherwise as unharmonized MEs.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the MEs stored in \code{seurat_obj@misc[[wgcna_name]]$hMEs} (harmonized) or \code{$MEs}.
#' @keywords scRNA-seq
#' @export
SetMEs <- function(seurat_obj, MEs, harmonized = TRUE, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # harmonized MEs?
  if (harmonized) {
    seurat_obj@misc[[wgcna_name]]$hMEs <- MEs
  } else {
    seurat_obj@misc[[wgcna_name]]$MEs <- MEs
  }
  seurat_obj
}


#' GetMEs
#'
#' Function to retrieve module eigengens from Seurat object.
#'
#' @param seurat_obj A Seurat object.
#' @param harmonized Logical; if TRUE (default), return harmonized MEs when available, otherwise unharmonized MEs.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame or matrix of module eigengenes (cells by modules). Falls back to unharmonized MEs if harmonized MEs are not present.
#' @keywords scRNA-seq
#' @export
GetMEs <- function(seurat_obj, harmonized = TRUE, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get harmonized MEs?
  if (harmonized == TRUE && !is.null(seurat_obj@misc[[wgcna_name]]$hMEs)) {
    MEs <- seurat_obj@misc[[wgcna_name]]$hMEs
  } else {
    MEs <- seurat_obj@misc[[wgcna_name]]$MEs
  }
  MEs
}


#' SetMELoadings
#'
#' @param seurat_obj A Seurat object.
#' @param loadings A named numeric vector of module eigengene feature loadings, appended to any existing loadings.
#' @param harmonized Logical; if TRUE (default), store as harmonized loadings, otherwise as unharmonized loadings.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the loadings appended to \code{seurat_obj@misc[[wgcna_name]]$hME_loadings} (harmonized) or \code{$ME_loadings}.
#' @keywords scRNA-seq
#' @export
SetMELoadings <- function(seurat_obj, loadings, harmonized = TRUE, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # harmonized MEs?
  if (harmonized) {
    seurat_obj@misc[[wgcna_name]]$hME_loadings <- c(seurat_obj@misc[[wgcna_name]]$hME_loadings, loadings)
  } else {
    seurat_obj@misc[[wgcna_name]]$ME_loadings <- c(seurat_obj@misc[[wgcna_name]]$ME_loadings, loadings)
  }
  seurat_obj
}


#' GetMELoadings
#'
#' Function to retrieve module eigengens from Seurat object.
#'
#' @param seurat_obj A Seurat object.
#' @param harmonized Logical; if TRUE (default), return harmonized loadings when available, otherwise unharmonized loadings.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A named vector of module eigengene feature loadings. Falls back to unharmonized loadings if harmonized loadings are not present.
#' @keywords scRNA-seq
#' @export
GetMELoadings <- function(seurat_obj, harmonized = TRUE, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get harmonized MEs?
  if (harmonized == TRUE && !is.null(seurat_obj@misc[[wgcna_name]]$hME_loadings)) {
    MEs <- seurat_obj@misc[[wgcna_name]]$hME_loadings
  } else {
    MEs <- seurat_obj@misc[[wgcna_name]]$ME_loadings
  }
  MEs
}


############################
# Module Scores
###########################


#' SetModuleScores
#'
#' @param seurat_obj A Seurat object.
#' @param mod_scores A data.frame of module expression scores (cells by modules).
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the module scores stored in \code{seurat_obj@misc[[wgcna_name]]$module_scores}.
#' @export
SetModuleScores <- function(seurat_obj, mod_scores, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$module_scores <- mod_scores
  seurat_obj
}


#' GetModuleScores
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of module expression scores (cells by modules), as computed by \code{ModuleExprScore}, or NULL if not computed.
#' @export
GetModuleScores <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$module_scores
}


############################
# Average Module Expression
###########################

#' SetAvgModuleExpr
#'
#' @param seurat_obj A Seurat object.
#' @param avg_mods A data.frame of the average expression of all genes in each module (cells by modules).
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the average module expression stored in \code{seurat_obj@misc[[wgcna_name]]$avg_modules}.
#' @export
SetAvgModuleExpr <- function(seurat_obj, avg_mods, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$avg_modules <- avg_mods
  seurat_obj
}


#' GetAvgModuleExpr
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of average module expression (cells by modules), or NULL if not computed.
#' @keywords scRNA-seq
#' @export
GetAvgModuleExpr <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$avg_modules
}


############################
# ModuleUMAP
###########################

#' SetModuleUMAP
#'
#' @param seurat_obj A Seurat object.
#' @param umap_df A data.frame of module UMAP coordinates for each gene.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the module UMAP stored in \code{seurat_obj@misc[[wgcna_name]]$module_umap}.
#' @keywords scRNA-seq
#' @export
SetModuleUMAP <- function(seurat_obj, umap_df, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$module_umap <- umap_df
  seurat_obj
}


#' GetModuleUMAP
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of module UMAP coordinates for each gene, as computed by \code{RunModuleUMAP}, or NULL if not computed.
#' @keywords scRNA-seq
#' @export
GetModuleUMAP <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$module_umap
}


############################
# ModuleTraitCorrelation
###########################


#' SetModuleTraitCorrelation
#'
#' @param seurat_obj A Seurat object.
#' @param mt_cor A list of module-trait correlation results with elements cor, pval, and fdr.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the module-trait correlations stored in \code{seurat_obj@misc[[wgcna_name]]$mt_cor}.
#' @keywords scRNA-seq
#' @export
SetModuleTraitCorrelation <- function(seurat_obj, mt_cor, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$mt_cor <- mt_cor
  seurat_obj
}


#' GetModuleTraitCorrelation
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A list with elements cor, pval, and fdr, each a list of module-by-trait matrices (one per group), as computed by \code{ModuleTraitCorrelation}.
#' @keywords scRNA-seq
#' @export
GetModuleTraitCorrelation <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$mt_cor
}


############################
# ModulePreservation
###########################


#' SetModulePreservation
#'
#' @param seurat_obj A Seurat object.
#' @param mod_pres The module preservation results to store (e.g. a list with elements obs and Z).
#' @param mod_name The name under which to store the module preservation results.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the results stored in \code{seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]]}.
#' @keywords scRNA-seq
#' @export
SetModulePreservation <- function(seurat_obj, mod_pres, mod_name, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # make an empty list if module preservation hasn't been called yet
  if (is.null(seurat_obj@misc[[wgcna_name]]$module_preservation)) {
    seurat_obj@misc[[wgcna_name]]$module_preservation <- list()
  }

  seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]] <- mod_pres
  seurat_obj
}


#' GetModulePreservation
#'
#' @param seurat_obj A Seurat object.
#' @param mod_name The name of the stored module preservation results to retrieve.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The module preservation results stored under \code{mod_name}; throws an error if none are found.
#' @keywords scRNA-seq
#' @export
GetModulePreservation <- function(seurat_obj, mod_name, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  if (is.null(seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]])) {
    stop('Invalid module preservation name.')
  }
  seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]]
}


############################
# GO term table
###########################

#' SetEnrichrTable
#'
#' @param seurat_obj A Seurat object.
#' @param enrich_table A data.frame of Enrichr results.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return The input Seurat object with the Enrichr table stored in \code{seurat_obj@misc[[wgcna_name]]$enrichr_table}.
#' @keywords scRNA-seq
#' @export
SetEnrichrTable <- function(seurat_obj, enrich_table, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set enrichr table
  seurat_obj@misc[[wgcna_name]]$enrichr_table <- enrich_table
  seurat_obj
}


#' GetEnrichrTable
#'
#' @param seurat_obj A Seurat object.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot. Defaults to the active experiment.
#' @return A data.frame of Enrichr results, as computed by \code{RunEnrichr}, or NULL if not computed.
#' @keywords scRNA-seq
#' @export
GetEnrichrTable <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$enrichr_table
}
