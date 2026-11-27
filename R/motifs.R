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
#' @export
MotifScan <- function(
  seurat_obj,
  pfm,
  EnsDb,
  species_genome,
  wgcna_name = NULL
) {
  check_suggests(
    c('motifmatchr', 'ensembldb', 'GenomeInfoDb', 'GenomicRanges', 'IRanges', 'BiocGenerics'),
    feature = 'MotifScan'
  )

  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get a dataframe of just the motif name and the motif ID:
  motif_df <- data.frame(
    motif_name = purrr::map(seq_along(pfm), function(i) {
      pfm[[i]]@name
    }) %>% unlist(),
    motif_ID = purrr::map(seq_along(pfm), function(i) {
      pfm[[i]]@ID
    }) %>% unlist()
  )

  # get promoter and gene coords:
  gene.promoters <- ensembldb::promoters(EnsDb)
  gene.coords <- ensembldb::genes(EnsDb)

  # subset by protein coding
  gene.promoters <- gene.promoters[gene.promoters$tx_biotype == 'protein_coding']
  gene.coords <- gene.coords[gene.coords$gene_biotype == 'protein_coding']

  # subset by main chromosomes
  gene.promoters <- gene.promoters[as.character(GenomeInfoDb::seqnames(gene.promoters)) %in% c(1:100, 'X', 'Y')]
  gene.coords <- gene.coords[as.character(GenomeInfoDb::seqnames(gene.coords)) %in% c(1:100, 'X', 'Y')]

  # add the gene name to the promoter object
  gene.promoters$symbol <- gene.coords$symbol[base::match(gene.promoters$gene_id, names(gene.coords))]

  # drop unnecessary chromosomes
  gene.promoters <- GenomeInfoDb::keepSeqlevels(
    gene.promoters,
    value = levels(droplevels(GenomeInfoDb::seqnames(gene.promoters)))
  )

  # rename seqlevels to add 'chr',
  old_levels <- levels(GenomeInfoDb::seqnames(gene.promoters))
  new_levels <- paste0('chr', old_levels)
  gene.promoters <- GenomeInfoDb::renameSeqlevels(gene.promoters, new_levels)

  # set the genome (not sure if we NEED to do this...)
  GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gene.promoters)) <- species_genome

  # set up promoters object that only has the necessary info for motifmatchr
  my_promoters <- GenomicRanges::GRanges(
    seqnames = droplevels(GenomeInfoDb::seqnames(gene.promoters)),
    IRanges::IRanges(
      start = BiocGenerics::start(gene.promoters),
      end = BiocGenerics::end(gene.promoters)
    ),
    symbol = gene.promoters$symbol,
    genome = species_genome
  )


  # scan these promoters for motifs:
  print('Matching motifs...')
  motif_ix <- motifmatchr::matchMotifs(pfm, my_promoters, genome = species_genome)

  # get the matches
  tf_match <- motifmatchr::motifMatches(motif_ix)
  rownames(tf_match) <- my_promoters$symbol

  # use motif names as the column names:
  # colnames(tf_match) <- motif_df$motif_name
  colnames(tf_match) <- motif_df$motif_ID

  #

  # only keep genes that are in the Seurat object and in the given EnsDb:
  gene_list <- rownames(seurat_obj)
  gene_list <- gene_list[gene_list %in% rownames(tf_match)]
  tf_match <- tf_match[gene_list, ]

  # get list of target genes for each TF:
  print('Getting putative TF target genes...')
  tf_targets <- list()
  n_targets <- list()
  for (cur_tf in colnames(tf_match)) {
    tf_targets[[cur_tf]] <- names(tf_match[, cur_tf][tf_match[, cur_tf]])
    n_targets[[cur_tf]] <- length(tf_targets[[cur_tf]])
  }
  n_targets <- unlist(n_targets)

  # add number of target genes to motif_df
  motif_df$n_targets <- n_targets

  # resolve cases where there's more than 1 motif per gene?
  # maybe take the union or smth.
  tfs_fix <- table(motif_df$motif_name)
  tfs_fix <- names(tfs_fix[which(tfs_fix > 1)])
  for (cur_tf in tfs_fix) {
    cur_motifs <- subset(motif_df, motif_name == cur_tf)
    cur_indices <- as.numeric(rownames(cur_motifs))
    cur_targets <- unique(unlist(tf_targets[cur_indices]))
    tf_targets[[cur_tf]] <- cur_targets
  }

  # add the gene name to the motif_df
  # remove extra characters from the motif names
  motif_names <- motif_df$motif_name
  tmp <- gsub('\\(.*)', '', motif_names)
  tmp <- gsub('::', ',', as.character(tmp))

  motif_df$tmp <- tmp

  # for motifs that correspond to two genes, split them apart
  tmp <- motif_df$tmp
  names(tmp) <- motif_df$motif_ID
  motif_df_tmp <- do.call(rbind, lapply(seq_along(tmp), function(i) {
    x <- tmp[i]
    id <- names(x)
    if (grepl(',', x)) {
      x <- as.character(unlist(do.call(rbind, strsplit(x, ','))))
    }
    data.frame(motif_ID = id, gene_name = as.character(x))
  }))

  # merge with the other motif df:
  ix <- match(motif_df_tmp$motif_ID, motif_df$motif_ID)
  motif_df_tmp <- cbind(motif_df_tmp, motif_df[ix, c('motif_name', 'n_targets')])
  rownames(motif_df_tmp) <- seq_len(nrow(motif_df_tmp))
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
OverlapModulesMotifs <- function(
  seurat_obj, wgcna_name = NULL
) {
  check_suggests('GeneOverlap', feature = 'OverlapModulesMotifs')

  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

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
  cur_module_genes <- modules %>%
    subset(module == cur_mod) %>%
    .$gene_name
  cur_targets <- as.character(unlist(tf_targets[1]))

  overlap_df <- do.call(rbind, lapply(mods, function(cur_mod) {
    cur_module_genes <- modules %>%
      subset(module == cur_mod) %>%
      .$gene_name
    cur_overlap_df <- do.call(rbind, lapply(tf_targets, function(cur_targets) {
      cur_overlap <- GeneOverlap::testGeneOverlap(GeneOverlap::newGeneOverlap(
        cur_module_genes,
        as.character(unlist(cur_targets)),
        genome.size = genome.size
      ))
      c(cur_overlap@odds.ratio, cur_overlap@pval, cur_overlap@Jaccard, length(cur_overlap@intersection))
    })) %>% as.data.frame()

    colnames(cur_overlap_df) <- c('odds_ratio', 'pval', 'Jaccard', 'size_intersection')
    cur_overlap_df$module <- cur_mod
    cur_overlap_df$tf <- names(tf_targets)

    # module color:
    cur_overlap_df$color <- modules %>%
      subset(module == cur_mod) %>%
      .$color %>%
      unique()
    cur_overlap_df
  }))

  # adjust for multiple comparisons:
  overlap_df$fdr <- p.adjust(overlap_df$pval, method = 'fdr')

  # significance level:
  overlap_df$Significance <- gtools::stars.pval(overlap_df$fdr)
  overlap_df$Significance <- ifelse(
    overlap_df$Significance == '.', '',
    overlap_df$Significance
  )

  # set factor levels for modules:
  overlap_df$module <- factor(overlap_df$module, levels = mods)

  # re-arrange columns:
  overlap_df <- overlap_df %>% dplyr::select(c(module, tf, color, odds_ratio, pval, fdr, Significance, Jaccard, size_intersection))

  # add overlap df to Seurat obj
  seurat_obj <- SetMotifOverlap(seurat_obj, overlap_df, wgcna_name)

  seurat_obj
}


#' Displays the top n TFs in a set of modules as a bar plot
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
MotifOverlapBarPlot <- function(
  seurat_obj,
  n_tfs = 10,
  plot_size = c(5, 6),
  outdir = 'MotifOverlaps/',
  motif_font = 'helvetica_regular',
  module_names = NULL,
  wgcna_name = NULL
) {
  check_suggests('ggseqlogo', feature = 'MotifOverlapBarPlot')

  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # make output dir if it doesn't exist:
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }

  # get Modules
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  if (is.null(module_names)) {
    module_names <- mods
  }

  # get overlap info from Seurat obj
  overlap_df <- GetMotifOverlap(seurat_obj, wgcna_name)
  motif_df <- GetMotifs(seurat_obj)

  # get pfm list from Seurat obj
  pfm <- GetPFMList(seurat_obj)

  # add motif ID to the overlap_df
  overlap_df$motif_ID <- motif_df$motif_ID[match(overlap_df$tf, motif_df$motif_name)]

  # subset by the modules that we are using:
  overlap_df <- overlap_df %>% subset(module %in% module_names)

  for (cur_mod in module_names) {
    print(cur_mod)
    # get df for cur mod
    plot_df <- overlap_df %>%
      subset(module == cur_mod) %>%
      top_n(n_tfs, wt = odds_ratio) %>%
      arrange(desc(odds_ratio))

    # make the barplot
    p1 <- plot_df %>% ggplot(aes(y = reorder(tf, odds_ratio), fill = odds_ratio, x = odds_ratio)) +
      geom_bar(stat = 'identity', width = 0.7) +
      NoLegend() +
      scale_fill_gradient(high = unique(plot_df$color), low = 'grey90') +
      ylab('') +
      theme(
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(
          t = 0, r = 0, b = 0, l = 0
        )
      )

    # make the motif logo plots:
    plot_list <- list()
    for (i in seq_len(nrow(plot_df))) {
      cur_id <- plot_df[i, 'motif_ID']
      cur_name <- plot_df[i, 'tf']
      plot_list[[cur_id]] <- ggplot() +
        ggseqlogo::geom_logo(as.matrix(pfm[[cur_id]]), font = motif_font) +
        ggseqlogo::theme_logo() +
        xlab('') +
        ylab(cur_name) +
        theme(
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(angle = 0),
          plot.margin = margin(
            t = 0, # Top margin
            r = 0, # Right margin
            b = 0, # Bottom margin
            l = 0
          ) # Left margin
        )
    }

    # wrap the motif logos
    patch1 <- wrap_plots(plot_list, ncol = 1)

    # assemble the plot with patchwork
    outplot <- (patch1 | p1) +
      plot_layout(ncol = 2, widths = c(1, 2)) +
      plot_annotation(
        title = paste0('Motif overlaps with ', cur_mod),
        theme = theme(plot.title = element_text(hjust = 0.5))
      )

    pdf(paste0(outdir, '/', cur_mod, '_motif_overlaps.pdf'), width = plot_size[1], height = plot_size[2], useDingbats = FALSE)
    print(outplot)
    dev.off()
  }
}
