
#' RunEnrichr
#'
#' Run Enrichr gene set enrichment tests on hdWGCNA modules
#'
#' @param seurat_obj A Seurat object
#' @param dbs character vector of EnrichR databases
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wait logical indicating whether or not to wait some time between sending requests to the EnrichR server.
#' @param wait_time the number of seconds to wait between sending requests to the EnrichR server. Value must be less than 60.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
#' RunEnrichr
RunEnrichr <- function(
  seurat_obj,
  dbs = c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021'),
  max_genes = 100,
  wait = TRUE,
  wait_time = 5,
  wgcna_name=NULL, ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # check wait_time 
  if(!is.numeric(wait_time)){
    stop(paste0('wait_time must be a numeric.'))
  } 
  if(wait_time > 60 | wait_time < 1){
    stop(paste0('Invalid value selected for wait_time, must be greater than 0 and less than 60.'))
  }

  # get modules:
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)

  # exclude grey
  mods <- mods[mods != 'grey']

  # hub genes 
  hub_df <- GetHubGenes(seurat_obj, n_hubs=max_genes, wgcna_name=wgcna_name)

  # run EnrichR for loop:
  combined_output <- data.frame()
  for(i in 1:length(mods)){
    
  	cur_mod <- mods[i] 
    cur_info <- subset(hub_df, module == cur_mod)
    cur_genes <- as.character(cur_info$gene_name)

    # run the enrichment test
    enriched <- enrichR::enrichr(cur_genes, dbs)

    if(wait){
      Sys.sleep(wait_time)
    }

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

#' RunEnrichrRegulons
#'
#' Run Enrichr gene set enrichment tests on hdWGCNA modules
#'
#' @param seurat_obj A Seurat object
#' @param dbs character vector of EnrichR databases
#' @param depth Include primary TF target genes (depth=1) or primary + secondary (depth=2)
#' @param use_regulons Use the regulons (default=TRUE) or the full set of TF target genes?
#' @param max_genes Max number of genes to include per module, ranked by kME.
#' @param wait logical indicating whether or not to wait some time between sending requests to the EnrichR server.
#' @param wait_time the number of seconds to wait between sending requests to the EnrichR server. Value must be less than 60.
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
RunEnrichrRegulons <- function(
    seurat_obj,
    dbs = c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021'),
    depth = 1,
    use_regulons=TRUE,
    min_genes = 5,
    wait = TRUE,
    wait_time = 5,
    wgcna_name=NULL, ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # check wait_time 
  if(!is.numeric(wait_time)){
    stop(paste0('wait_time must be a numeric.'))
  } 
  if(wait_time > 60 | wait_time < 1){
    stop(paste0('Invalid value selected for wait_time, must be greater than 0 and less than 60.'))
  }

  # get regulons
  tf_regulons <- GetTFRegulons(seurat_obj, wgcna_name)
  tf_list <- unique(tf_regulons$tf)

  # run EnrichR for loop:
  combined_output <- data.frame()
  for(i in 1:length(tf_list)){
  	cur_tf <- tf_list[i]
    cur_targets <- GetTFTargetGenes(seurat_obj, cur_tf, depth=depth, use_regulons=use_regulons, wgcna_name=wgcna_name)
    cur_pos <- subset(cur_targets, Cor > 0) %>% .$gene %>% unique()
    cur_neg <- subset(cur_targets, Cor < 0) %>% .$gene %>% unique()

    # run the enrichment test
    if(length(cur_pos) >= min_genes){
        enriched <- enrichR::enrichr(cur_pos, dbs)

        # collapse into one db
        for(db in names(enriched)){
            cur_df <- enriched[[db]]
            if (nrow(cur_df) > 1){
                cur_df$db <- db
                cur_df$tf <- cur_tf
                cur_df$target_type <- 'positive'
                combined_output <- rbind(combined_output, cur_df)
            }
        }

        if(wait){
            Sys.sleep(wait_time)
        }
    }

    if(length(cur_neg) >= min_genes){
        enriched <- enrichR::enrichr(cur_neg, dbs)

        # collapse into one db
        for(db in names(enriched)){
            cur_df <- enriched[[db]]
            if (nrow(cur_df) > 1){
                cur_df$db <- db
                cur_df$tf <- cur_tf
                cur_df$target_type <- 'negative'
                combined_output <- rbind(combined_output, cur_df)
            }
        }

        if(wait){
            Sys.sleep(wait_time)
        }
    }

  }

  # add GO term table to seurat object
  seurat_obj <- SetEnrichrRegulonTable(seurat_obj, combined_output, wgcna_name)
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

