
FindMajorIsoforms <- function(
  seurat_obj,
  group.by,
  replicate_col,
  isoform_delim = '[.]',
  proportion_thresh = 0.8,
  low_thresh = 25,
  assay = 'iso',
  slot = 'counts',
  wgcna_name = NULL
){

  # set as active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

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

  # TODO: add checks for this step?
  # create dataframe to map isoforms to gene names
  if(is.null(iso_df)){
    iso_names <- rownames(seurat_obj)
    gene_names <- do.call(rbind, strsplit(iso_names, iso_delim))[,1]
    iso_df <- data.frame(
      iso = iso_names,
      gene = gene_names
    )
  }
  if(max(table(iso_df$iso)) > 1){
    stop('Each isoform must only appear once in iso_df. There is likely an issue with the isoform_delim.')
  }

  # get counts matrix:
  X <- Seurat::GetAssayData(seurat_obj, assay=assay, slot=slot)

  # get pseudobulk replicates
  pseudobulks <- Libra::to_pseudobulk(
    X, meta = seurat_obj@meta.data,
    cell_type_col = group.by,
    replicate_col = replicate_col,
    label_col = replicate_col,
    min_reps=0
  )
  
  # loop through each group to find the major isoforms
  major_list <- list()
  for(cur_group in names(pseudobulks)){
    print(cur_group)
    cur_pb <- pseudobulks[[cur_group]]
    cur_pb_genes <- do.call(rbind, strsplit(rownames(cur_pb), iso_delim))[,1]

    # loop through each gene:
    major_isos <- sapply(unique(iso_df$gene), function(cur_gene){
      cur_iso_df <- subset(iso_df, gene == cur_gene)
      if(nrow(cur_iso_df) == 1){
        return(cur_iso_df$iso)
      }

      cur_iso_pb <- cur_pb[cur_pb_genes == cur_gene,]
      cur_iso_sum  <- rowSums(cur_iso_pb)
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

  }

  major_list

}