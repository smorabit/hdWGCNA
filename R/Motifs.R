#' MotifScan
#'
#' @description
#' This function scans the promoter regions of protein-coding genes for transcription factor (TF) motifs.
#' It extracts promoter sequences using an Ensembl database (`EnsDb`) and then searches these sequences
#' for TF binding motifs using position frequency matrices (PFMs). The Seurat object is updated with a matrix
#' of motif-gene matches, a list of target genes for each TF, and additional motif information.
#'
#' @param seurat_obj A Seurat object that will be updated with motif scan results.
#' @param pfm A list of position frequency matrices (PFMs), such as those from the JASPAR2020 database.
#' @param EnsDb An Ensembl database object (e.g., `EnsDb.Hsapiens.v86` or `EnsDb.Mmusculus.v79`) containing gene and promoter annotations.
#' @param species_genome A character string specifying the genome version (e.g., "hg38" for human or "mm10" for mouse).
#' @param wgcna_name A character string specifying the name of the WGCNA experiment to associate with the motif data (optional).
#' If NULL, the active WGCNA experiment in `seurat_obj@misc` will be used.
#'
#' @details 
#' The `MotifScan` function performs the following steps:
#' - It extracts promoter regions (typically 2 kb upstream of the transcription start site) of protein-coding genes from the provided `EnsDb`.
#' - It uses the `motifmatchr` package to search these promoters for TF motifs, using the input PFMs.
#' - The function returns the Seurat object updated with several key pieces of information:
#'   - A motif-gene match matrix indicating the presence or absence of each motif in the promoter of each gene.
#'   - A list of target genes for each TF based on motif presence.
#'   - A summary of motifs, including the number of target genes for each TF.
#' - This data is stored in the `seurat_obj`'s metadata and can be used for downstream analysis, such as regulatory network inference.
#'
#' @return A modified Seurat object containing the results of the motif scan, including:
#' - `seurat_obj@misc$motif_matrix`: A binary matrix indicating motif matches for each gene.
#' - `seurat_obj@misc$motif_info`: A data frame containing motif names, IDs, and the number of target genes.
#' - `seurat_obj@misc$motif_targets`: A list of genes targeted by each motif.
#' - `seurat_obj@misc$pfm`: The original PFMs used for the motif scan.
#'
#' @import Seurat 
#' @importFrom ensembldb promoters genes
#' @importFrom GenomeInfoDb seqnames seqlevels keepSeqlevels renameSeqlevels
#' @importFrom GRanges GRanges
#' @importFrom motifmatchr matchMotifs motifMatches
#' @export
MotifScan <- function(
    seurat_obj,
    pfm, 
    EnsDb, 
    species_genome, 
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # get a dataframe of just the motif name and the motif ID:
    motif_df <- data.frame(
        motif_name = purrr::map(1:length(pfm), function(i){pfm[[i]]@name}) %>% unlist,
        motif_ID = purrr::map(1:length(pfm), function(i){pfm[[i]]@ID}) %>% unlist
    )

    # get promoter and gene coords:
    gene.promoters <- ensembldb::promoters(EnsDb)
    gene.coords <- ensembldb::genes(EnsDb)

    # subset by protein coding
    gene.promoters <- gene.promoters[gene.promoters$tx_biotype == 'protein_coding']
    gene.coords <- gene.coords[gene.coords$gene_biotype == 'protein_coding']

    # subset by main chromosomes
    gene.promoters <- gene.promoters[as.character(GenomeInfoDb::seqnames(gene.promoters)) %in% c(1:100, 'X','Y')]
    gene.coords <- gene.coords[as.character(GenomeInfoDb::seqnames(gene.coords)) %in% c(1:100, 'X', 'Y')]

    # add the gene name to the promoter object
    gene.promoters$symbol <- gene.coords$symbol[base::match(gene.promoters$gene_id, names(gene.coords))]

    # drop unnecessary chromosomes
    gene.promoters <- GenomeInfoDb::keepSeqlevels(
        gene.promoters, 
        value=levels(droplevels(GenomeInfoDb::seqnames(gene.promoters)))
    )

    # rename seqlevels to add 'chr', 
    old_levels <- levels(GenomeInfoDb::seqnames(gene.promoters))
    new_levels <- paste0('chr', old_levels)
    gene.promoters <- GenomeInfoDb::renameSeqlevels(gene.promoters, new_levels)

    # set the genome (not sure if we NEED to do this...)
    GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gene.promoters)) <- species_genome

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

    # add number of target genes to motif_df
    motif_df$n_targets <- n_targets

    # add the gene name to the moti_df
    # remove extra characters from the motif names
    motif_names <- motif_df$motif_name
    tmp <- gsub("\\(.*)", "", motif_names)
    tmp <- gsub('::', ',', as.character(tmp))

    motif_df$tmp <- tmp

    # for motifs that correspond to two genes, split them apart
    tmp <- motif_df$tmp; names(tmp) <- motif_df$motif_ID
    motif_df_tmp <- do.call(rbind,lapply(1:length(tmp), function(i){
    x <- tmp[i]
    id <- names(x)
    if(grepl(',', x)){
        x <- as.character(unlist(do.call(rbind, strsplit(x, ','))))
    }
    data.frame(motif_ID = id, gene_name = as.character(x))
    }))

    # merge with the other motif df:
    ix <- match(motif_df_tmp$motif_ID, motif_df$motif_ID)
    motif_df_tmp <- cbind(motif_df_tmp, motif_df[ix,c('motif_name', 'n_targets')])
    rownames(motif_df_tmp) <- 1:nrow(motif_df_tmp)
    motif_df_tmp <- dplyr::select(motif_df_tmp, c(motif_ID, motif_name, n_targets, gene_name))

    motif_df <- motif_df_tmp

    # subset to only contain genes in the seurat obj
    motif_df <- subset(motif_df, gene_name %in% rownames(seurat_obj))

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


