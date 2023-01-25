
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
