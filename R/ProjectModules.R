
#' ProjectModules
#'
#' Project a set of co-expression modules from a reference to a query dataset
#'
#' @param seurat_obj A Seurat object where the modules will be projected to
#' @param seurat_ref A Seurat object containing the co-expression modules to be projected
#' @param modules optionally provide a dataframe containing gene-module information to bypass seurat_ref
#' @param group.by.vars groups to harmonize by
#' @param gene_mapping dataframe to map gene names between genomes
#' @param genome1_col column in gene_mapping containing the genes for seurat_ref (reference)
#' @param genome2_col column in gene_mapping containing the genes for seurat_obj (query)
#' @param overlap_proportion the proportion of genes that must be present in seurat_obj for a given module to be projected. Default = 0.5 (50% of genes)
#' @param vars.to.regress character vector of variables in seurat_obj@meta.data to regress when running ScaleData
#' @param scale.model.use model to scale data when running ScaleData choices are "linear", "poisson", or "negbinom"
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_ref@misc slot
#' @param wgcna_name_proj The name of the hdWGCNA experiment to be created for the projected modules in seurat_obj
#' @keywords scRNA-seq
#' @export
#' @examples
#' ProjectModules
ProjectModules <- function(
  seurat_obj,
  seurat_ref,
  modules=NULL,
  group.by.vars=NULL,
  gene_mapping=NULL, # table mapping genes from species 1 to species 2
  genome1_col=NULL, genome2_col=NULL,
  overlap_proportion = 0.5,
  vars.to.regress = NULL,
  scale.model.use = 'linear',
  wgcna_name=NULL, wgcna_name_proj=NULL,
  ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_ref@misc$active_wgcna}
  CheckWGCNAName(seurat_ref, wgcna_name)

  # get modules to be projected:
  if(is.null(modules)){
    modules <- GetModules(seurat_ref, wgcna_name)
  } else{
    if(!all(c("gene_name", "module", "color") %in% colnames(modules))){
      stop('Missing columns in modules table. Required columns are gene_name, module, color')
    }
  }

  # cast "modules" to a factor
  if(!is.factor(modules$module)){
    modules$module <- as.factor(modules$module)
  }

  # are we mapping gene names?
  if(!is.null(gene_mapping)){
    modules <- TransferModuleGenome(modules, gene_mapping, genome1_col, genome2_col)
  }

  # what is the proportion of genes in each module that are in seurat_obj?
  mods <- as.character(unique(modules$module))
  mod_props <- unlist(lapply(mods, function(x){
    cur_mod <- subset(modules, module == x)
    sum(cur_mod$gene_name %in% rownames(seurat_obj)) / nrow(cur_mod)
  }))
  mods_keep <- mods[mod_props >= overlap_proportion]
  mods_remove <- mods[mod_props < overlap_proportion]

  if(length(mods) > 0){
    warning(paste0("The following modules will not be projected because too few genes are present in seurat_obj: ", paste(mods_remove, collapse=', ')))
  }

  # only keep modules that have enough overlapping genes
  modules <- subset(modules, module %in% mods_keep) %>%
    dplyr::mutate(module = droplevels(module))

  # get genes that overlap between WGCNA genes & seurat_obj genes:
  gene_names <- modules$gene_name
  genes_use <- intersect(gene_names, rownames(seurat_obj))

  # subset modules by genes in this seurat object:
  modules <- modules %>% subset(gene_name %in% genes_use)

  # setup new seurat obj for wgcna:
  if(!(wgcna_name_proj %in% names(seurat_obj@misc))){
    seurat_obj <- SetupForWGCNA(
      seurat_obj,
      wgcna_name = wgcna_name_proj,
      features = genes_use
    )
  }

  # project modules:
  print('project modules')
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars=group.by.vars,
    modules = modules,
    vars.to.regress = vars.to.regress,
    scale.model.use = scale.model.use,
    wgcna_name = wgcna_name_proj,
    ...
  )

  seurat_obj

}



#' TransferModuleGenome
#'
#' Takes a module table and a gene mapping table (like from biomart) with gene
#' names from two genomes in order to switch
#'
#' @param modules module table from GetModules function
#' @param gene_mapping table that has the gene names for both genomes
#' @param genome1_col the column in gene_mapping that has gene names for the genome currently present in modules
#' @param genome2_col the column in gene_mapping that has gene names for the genome to transfer to
#' @keywords scRNA-seq
#' @export
#' @examples
#' TransferModuleGenome
TransferModuleGenome <- function(
  modules, gene_mapping,
  genome1_col=NULL, genome2_col=NULL
){

  # use the first & second columns if these are null
  if(is.null(genome1_col)){genome1_col <- colnames(gene_mapping)[1]}
  if(is.null(genome2_col)){genome2_col <- colnames(gene_mapping)[2]}

  # need to add a check that there's a one to one mapping
  if(!all(table(gene_mapping[[genome2_col]]) == 1)){
    stop("There can only be one value for each feature in genome2_col in the gene_mapping table.")
  }

  # switch gene names to human:
  gene_mapping <- gene_mapping[,c(genome1_col, genome2_col)]

  # only keep genome1 genes that are in the WGCNA gene list:
  gene_list <- modules$gene_name
  gene_mapping <- gene_mapping[gene_mapping[,genome1_col] %in% gene_list,]

  # subset modules
  modules <- subset(modules, gene_name %in% gene_mapping[,genome1_col])

  # match the order of the given gene list, and remove NA entries
  gene_match <- match(gene_list, gene_mapping[,genome1_col])
  gene_mapping <- na.omit(gene_mapping[gene_match,])

  # update modules table with the new gene names
  # modules <- na.omit(modules[gene_match,])

  print(head(gene_mapping))

  modules$gene_name <- gene_mapping[,genome2_col]
  rownames(modules) <- modules$gene_name

  modules

}
