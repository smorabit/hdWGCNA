#' RunModuleUMAP
#'
#' Run UMAP on co-expression matrix using hub genes as features.
#'
#' @param seurat_obj A Seurat object
#' @param n_hubs number of hub genes to use in the UMAP computation
#' @param exclude_grey logical indicating whether to include grey module
#' @param genes_use character vector of genes to use for the UMAP, must already be present in GetModules(seurat_obj)
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @param ... Additional parameters supplied to uwot::umap
#' @import uwot
#' @import Seurat
#' @export
RunModuleUMAP <- function(
  seurat_obj,
  n_hubs = 10,
  exclude_grey = TRUE,
  genes_use = NULL,
  wgcna_name = NULL,
  n_neighbors= 25,
  metric = "cosine",
  spread=1,
  min_dist=0.4,
  supervised = FALSE,
  ...
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get modules,
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)

  # check if we have eigengene-based connectivities:
  if(!all(paste0('kME_', as.character(mods)) %in% colnames(modules))){
    stop('Eigengene-based connectivity (kME) not found. Did you run ModuleEigengenes and ModuleConnectivity?')
  }

  # handle excluding grey module
  if(exclude_grey){
    mods <- mods[mods != 'grey']
    modules <- subset(modules, module != 'grey')
  }

  # get hub genes:
  hub_df <- GetHubGenes(seurat_obj, n_hubs=n_hubs, wgcna_name=wgcna_name)

  # get all genes that aren't in gray mod
  selected_genes <- modules[modules$module %in% mods,'gene_name']

  if(!is.null(genes_use)){
    selected_genes <- selected_genes[selected_genes %in% genes_use]
  } 

  # get the TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # subset the TOM for umap
  # keep all genes as rows, and keep only hubs as cols (features)
  feature_mat <- TOM[selected_genes,hub_df$gene_name]

  # run UMAP
  if(supervised){
    hub_umap <-  uwot::umap(
      X = feature_mat,
      min_dist = min_dist,
      n_neighbors= n_neighbors,
      metric = metric,
      spread=spread,
      y = modules$module, # for supervised UMAP
      ...
    )
  } else {
    hub_umap <-  uwot::umap(
      X = feature_mat,
      min_dist = min_dist,
      n_neighbors= n_neighbors,
      metric = metric,
      spread=spread,
      ...
    )
  }

  # set up plotting df
  plot_df <- as.data.frame(hub_umap)
  colnames(plot_df) <- c("UMAP1", "UMAP2")
  plot_df$gene <- rownames(feature_mat)

  # add module color, and hub gene status to the plotting df:
  ix <- match(plot_df$gene, modules$gene_name)
  plot_df$module <- modules$module[ix]
  plot_df$color <- modules$color[ix]
  plot_df$hub <- ifelse(
    plot_df$gene %in% as.character(hub_df$gene_name), 'hub', 'other'
  )

  # get kME values for each gene
  kMEs <- do.call(rbind, lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur <- cur[,c('gene_name', paste0('kME_', cur_mod))]
    colnames(cur) <- c('gene_name', 'kME')

    # scale kMEs between 0 & 1:
    cur$kME <- scale01(cur$kME)
    cur
  }))
  ix <- kMEs$gene_name[match(plot_df$gene, kMEs$gene_name)]
  plot_df$kME <- kMEs[ix, 'kME']

  # add the UMAP to the Seurat object:
  seurat_obj <- SetModuleUMAP(seurat_obj, plot_df, wgcna_name)

  seurat_obj
}
