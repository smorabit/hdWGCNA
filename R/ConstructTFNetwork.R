#' ConstructTFNetwork
#'
#' Construct a network of transcription factors and target genes based on gene co-expression
#'
#' @return seurat_obj with the TF network added
#'
#' @param seurat_obj A Seurat object
#' @param model_params a list of model parameters to pass to xgboost
#' @param nfold number of folds for cross validation
#' @param wgcna_name name of the WGCNA experiment
#' 
#' @details 
#' ConstructTFNetwork uses motif-gene information to build a directed network of transcription 
#' factors (TFs) and target genes. XGBoost regression is leveraged to model the expression of 
#' each gene based on its candidate TF regulators. This analysis gives us information about 
#' how important each TF is for predicting each gene, allowing us to prioritize the most likely
#' regulators of each gene. This process is done on the gene expression matrix stored with SetDatExpr,
#' which is typically the hdWGCNA metacell gene expression matrix.  
#' 
#' @import Seurat, Matrix, xgboost
#' @export
ConstructTFNetwork <- function(
    seurat_obj,
    model_params,
    nfold=5,
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # get the motif information from the Seurat object:
    motif_matrix <- GetMotifMatrix(seurat_obj)
    motif_df <- GetMotifs(seurat_obj)

    if(is.null(motif_df)){
        stop("Motif info not found in the seruat_obj. Please run MotifScan first.")
    }

    # check that gene names column is in the motif names
    if(! 'gene_name' %in% colnames(motif_df)){
        stop('gene_name column missing in motif table (GetMotifs(seurat_obj)). Please add a column indicating the gene_name in the seurat_obj for each motif.' )
    }

    # subset the motif_df by genes that are in the seurat obj:
    motif_df <- subset(motif_df, gene_name %in% rownames(seurat_obj))

    # get the expression matrix:
    datExpr <- as.matrix(GetDatExpr(seurat_obj, wgcna_name=wgcna_name))
    genes_use <- colnames(datExpr)

    # set up output dataframes
    importance_df <- data.frame()
    eval_df <- data.frame()

    # set up the progress bar
    pb <- utils::txtProgressBar(min = 0, max = length(genes_use), style = 3, width = 50, char = "=")
    counter <- 1

    for(cur_gene in genes_use){

        setTxtProgressBar(pb, counter)

        # check if this gene is in the motif matrix:
        if(! cur_gene %in% rownames(motif_matrix)){
            print(paste0('Gene not found in the motif_matrix, skipping ', cur_gene))
            next
        }

        # get the list of TFs that regulate this gene:
        cur_tfs <- names(which(motif_matrix[cur_gene,]))
        cur_tfs <- subset(motif_df, motif_name %in% cur_tfs) %>% .$gene_name %>% unique
        cur_tfs <- cur_tfs[cur_tfs %in% genes_use]

        # set up the expression matrices
        if(cur_gene %in% cur_tfs){
            cur_tfs <- cur_tfs[cur_tfs != cur_gene]
            x_vars <- datExpr[,cur_tfs]
        } 
        x_vars <- datExpr[,cur_tfs]
        y_var <- as.numeric(datExpr[,cur_gene])

        if(length(cur_tfs) < 2){
            print(paste0('Not enough putative TFs, skipping ', cur_gene))
            next
        }

        # correlation:
        tf_cor <- as.numeric(cor(
            x=as.matrix(x_vars),
            y=y_var
        ))
        names(tf_cor) <- cur_tfs

        # run xgboost model
        if(all(y_var == 0)){
            print(paste0('skipping ', cur_gene))
            next
        }
        xgb <- xgboost::xgb.cv(
            params = model_params,
            data = x_vars,
            label = y_var,
            nrounds = 100,
            showsd = FALSE,
            nfold = nfold,
            callbacks = list(cb.cv.predict(save_models=TRUE)),
            verbose=FALSE
        )

        # get the CV evaluation info
        xgb_eval <- as.data.frame(xgb$evaluation_log)
        xgb_eval$variable <- cur_gene

        # average the importance score from each fold
        importance <- Reduce('+', lapply(1:nfold, function(i){
            cur_imp <- xgb.importance(feature_names = colnames(x_vars), model = xgb$models[[i]])
            ix <- match(colnames(x_vars),  as.character(cur_imp$Feature))
            cur_imp <- as.matrix(cur_imp[ix,-1])
            cur_imp[is.na(cur_imp)] <- 0
            cur_imp
        })) / nfold
        importance <- as.data.frame(importance)

        # add tf and source info
        importance$tf <- colnames(x_vars)
        importance$gene<- cur_gene

        # add the tf correlation information
        importance$Cor <- as.numeric(tf_cor)

        # re-order columns, and re-order rows by gain:
        importance <- importance %>% dplyr::select(c(tf, gene, Gain, Cover, Frequency, Cor))
        importance <- arrange(importance, -Gain)

        # append
        importance_df <- rbind(importance_df, importance)
        eval_df <- rbind(eval_df, xgb_eval)

        # update progress bar
        counter <- counter+1

    } 

    # close the progress bar
    close(pb)

    # add the results to the Seurat object
    seurat_obj <- SetTFNetwork(seurat_obj, importance_df, wgcna_name=wgcna_name)
    seurat_obj <- SetTFEval(seurat_obj, eval_df, wgcna_name=wgcna_name)

    #return(list(importance=importance_df, eval=eval_df))
    seurat_obj

}