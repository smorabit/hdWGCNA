#' FindMajorIsoforms
#'
#' Finds the set of "major isoforms" for each gene in each cell population. The set of 
#' major isoforms consists of the isoforms accounting for a desired proportion of a gene's 
#' expression. Pseudobulk replicates in each sample and cell population are used for this calculation.
#' 
#' @return a list of major isoforms for each cell population
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing grouping info, ie clusters or celltypes
#' @param replicate_col column in seurat_obj@meta.data denoting each replicate / sample
#' @param isoform_delim 
#' @param proportion_thresh desired proportion of expression to define the set of major isoforms. Default = 0.8.
#' @param low_thresh lower bound of expression level for considering an isoform as part of the major isoform set. 
#' @param assay Assay in seurat_obj containing isoform expression information.
#' @param slot Slot in seurat_obj, default to counts slot.
#' @param cluster_markers Cell population marker gene table from Seurat FindAllMarkers for the same cell populations specified in group.by. Optional parameter, will exclude isoforms that are not from marker genes.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @details
#' FindMajorIsoforms computes the set of major isoforms in a given Seurat object that contains isoform-level 
#' expression information. First, pseudobulk replicates are computed for the given cell populations and samples 
#' present in the Seurat object. For each gene in each cell population, we rank the gene's isoforms by expression 
#' level and take the top expressing isoforms that make up the desired proportion of the gene's total expression,
#' making sure to exclude any very lowly expressed isoforms.
#' 
#' Optionally, the user may supply a marker gene table for each cell population (formatted like the output of Seurat 
#' FindAllMarkers), and then the algorithm will only return major isoforms of the given marker genes.
#' @import Seurat
#' @import Matrix
#' @import magrittr
#' @export
FindMajorIsoforms <- function(
  seurat_obj,
  group.by,
  replicate_col,
  isoform_delim = '[.]',
  proportion_thresh = 0.8,
  low_thresh = 25,
  assay = 'iso',
  slot = 'counts',
  cluster_markers = NULL,
  wgcna_name = NULL
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

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

  # check that replicate_col is valid
  if(!(class(seurat_obj@meta.data[,replicate_col]) %in% c('character', 'factor'))){
    stop('Selected replicate_col must be a character or a factor, but ', replicate_col, ' is a ', class(seurat_obj@meta.data[,replicate_col]), '.')
  }

  # check that the current assay is the selected assay:
  if(DefaultAssay(seurat_obj) != assay){
    stop('DefaultAssay(seurat_obj) is not ', assay, ' please switch the default assay to the desired assay before running this function.')
  }

  # check the marker gene table:
  if(!is.null(cluster_markers)){
    # check columns 
    if(!all(c("p_val", 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster', 'gene') %in% colnames(cluster_markers))){
      stop('Invalid format for cluster_markers table. Expecting a table from Seurat FindAllMarkers, must have the following columns present in the table: p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene.')
    }

    # check the clusters:
    if(!all(as.character(unique(seurat_obj@meta.data[,group.by])) %in% as.character(unique(cluster_markers$cluster)))){
      stop('Groups in group.by must also be present in cluster_markers$cluster')
    }
  }

  # create dataframe to map isoforms to gene names
  iso_names <- rownames(seurat_obj)
  gene_names <- do.call(rbind, strsplit(iso_names, isoform_delim))[,1]
  iso_df <- data.frame(
    iso = iso_names,
    gene = gene_names
  )
  if(max(table(iso_df$iso)) > 1){
    stop('Each isoform must only appear once in iso_df. There is likely an issue with the isoform_delim.')
  }

  # get counts matrix:
  X <- Seurat::GetAssayData(seurat_obj, assay=assay, slot=slot)
 
  # get pseudobulk replicates
  message("Constructing pseudobulks...")
  pseudobulks <- to_pseudobulk(
    X, meta = seurat_obj@meta.data,
    cell_type_col = group.by,
    replicate_col = replicate_col,
    label_col = replicate_col,
    min_reps=0
  )
  message("Done!")

  # initialize progress bar
  message("Finding major isoforms...")
  pb <- utils::txtProgressBar(min = 0, max = length(pseudobulks) ,style = 3, width = 50, char = "=")
  pb_counter <- 1

  # loop through each group to find the major isoforms
  major_list <- list()
  for(cur_group in names(pseudobulks)){
    #print(cur_group)
    cur_pb <- pseudobulks[[cur_group]]
    cur_pb_genes <- do.call(rbind, strsplit(rownames(cur_pb), isoform_delim))[,1]

    # loop through each gene:
    major_isos <- sapply(unique(iso_df$gene), function(cur_gene){
      cur_iso_df <- subset(iso_df, gene == cur_gene)
      if(nrow(cur_iso_df) == 1){
        return(cur_iso_df$iso)
      }

      # get the pseudobulk expression for all the isoforms in this gene
      cur_iso_pb <- cur_pb[cur_pb_genes == cur_gene,]
      cur_iso_sum  <- rowSums(cur_iso_pb)

      # compute the proportion of the gene's total expression in each isoform
      cur_iso_prop <- cur_iso_sum / sum(cur_iso_sum)

      # don't consider isoforms with very few counts
      ix <- cur_iso_sum > low_thresh
      cur_iso_sum <- cur_iso_sum[ix]
      if(length(cur_iso_sum) == 0){
        return()
      }

      # order isoforms from high to low:
      cur_iso_prop <- cur_iso_prop[rev(order(cur_iso_prop))]

      # identify the set of genes that reach the chosen propotion
      prop_sums <- sapply(1:length(cur_iso_prop), function(i){
        sum(cur_iso_prop[1:i])
      })
      ix <- min(which(prop_sums >= proportion_thresh))

      # return the set of major isoforms
      names(cur_iso_prop[cur_iso_prop >= cur_iso_prop[ix]])

    })
    major_list[[cur_group]] <- as.character(unlist(major_isos))

    # update progress bar
    setTxtProgressBar(pb, pb_counter)
    pb_counter <- pb_counter+1

  }

  # close progress bar
  close(pb)
  message("Done!")

  # return here if the cluster marker table isn't specified
  if(is.null(cluster_markers)){
    return(major_list)
  }

  # intersect with the markers table if it's provided:
  major_marker_list <- lapply(names(major_list), function(cur_group){
    cur_isos <- major_list[[cur_group]]
    cur_genes <- unique(do.call(rbind, strsplit(cur_isos, isoform_delim))[,1])
    cur_genes <- subset(cluster_markers, cluster == cur_group & gene %in% cur_genes) %>% .$gene
    subset(iso_df, gene %in% cur_genes & iso %in% cur_isos) %>% .$iso
  })
  names(major_marker_list) <- names(major_list)

  # return major isoforms intersected with the marker gene table
  list('major' = major_list, 'marker' = major_marker_list)

}



