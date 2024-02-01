

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

  # TODO: add checks?
  # check pfm
  # check ensdb
  
  print('set up motif df')

  # get a dataframe of just the motif name and the motif ID:
  motif_df <- data.frame(
    motif_name = purrr::map(1:length(pfm), function(i){pfm[[i]]@name}) %>% unlist,
    motif_ID = purrr::map(1:length(pfm), function(i){pfm[[i]]@ID}) %>% unlist
  )

  print('get promoters from ensdb')

  # get promoters of protein coding genes from the given Ensdb
  # note: everything breaks if I try to use X & Y chromosomes.
 
  # get promoter and gene coords:
  gene.promoters <- ensembldb::promoters(EnsDb)
  gene.coords <- ensembldb::genes(EnsDb)
  print('here')

  # subset by protein coding
  gene.promoters <- gene.promoters[gene.promoters$tx_biotype == 'protein_coding']
  gene.coords <- gene.coords[gene.coords$gene_biotype == 'protein_coding']
  print('here')

  # subset by chromosomes
  gene.promoters <- gene.promoters[as.character(GenomeInfoDb::seqnames(gene.promoters)) %in% c(1:100, 'X','Y')]
  gene.coords <- gene.coords[as.character(GenomeInfoDb::seqnames(gene.coords)) %in% c(1:100, 'X', 'Y')]
  print('here')


  # gene.promoters <- ensembldb::promoters(EnsDb, filter = ~ gene_biotype == "protein_coding") 
  # gene.promoters <- gene.promoters[GenomeInfoDb::seqnames(gene.promoters) %in% c(1:100)]
  
  
  # gene.coords <- ensembldb::genes(EnsDb)

  # gene.coords <- ensembldb::genes(EnsDb, filter = ~ gene_biotype == "protein_coding")
  # print('here')

  # gene.coords <- gene.coords[GenomeInfoDb::seqnames(gene.coords) %in% c(1:100)]

  print('add gene name to promoter')

  # add the gene name to the promoter object
  gene.promoters$symbol <- gene.coords$symbol[base::match(gene.promoters$gene_id, names(gene.coords))]

  print('drop chromosomes')

  # drop unnecessary chromosomes
  gene.promoters <- GenomeInfoDb::keepSeqlevels(
    gene.promoters, 
    value=levels(droplevels(GenomeInfoDb::seqnames(gene.promoters)))
  )

  print('rename seqlevels')

  # rename seqlevels to add 'chr', 
  old_levels <- levels(GenomeInfoDb::seqnames(gene.promoters))
  #new_levels <- ifelse(old_levels %in% c('X', 'Y'), old_levels, paste0('chr', old_levels))
  new_levels <- paste0('chr', old_levels)
  gene.promoters <- GenomeInfoDb::renameSeqlevels(gene.promoters, new_levels)

  print('set the genome')

  # set the genome (not sure if we NEED to do this...)
  GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gene.promoters)) <- species_genome

  print('set up the promoters granges object')

  # set up promoters object that only has the necessary info for motifmatchr
  my_promoters <- GRanges(
    seqnames =  droplevels(GenomeInfoDb::seqnames(gene.promoters)),
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
  print('Getting putative TF target genes...')
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
