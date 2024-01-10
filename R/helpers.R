

#' scale01
#'
#' Function to scale a numeric vector between 0 and 1.
#'
#' @param x numeric vector
#' @keywords scRNA-seq
#' @export
#' @examples scale01
scale01 <- function(x){
  y <- min(x); z <- max(x);
  (x-y) / (z-y)
}



#' umap_theme
#'
#' ggplot theme to remove axes etc.
#'
#' @keywords scRNA-seq
#' @export
#' @examples umap_theme
umap_theme <- function(){
  theme(
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
}



#' shuffle_points
#'
#' Function to shuffle the rows of a dataframe.
#'
#' @param df dataframe
#' @keywords scRNA-seq
#' @export
#' @examples shuffle_points
shuffle_points <- function(df){
  return(df[sample(1:dim(df)[1], dim(df)[1]),])
}


#' Sparse matrix correlation
#' 
#' Compute the Pearson correlation matrix between
#' columns of two sparse matrices.
#' 
#' Originally from 
#' \url{http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r}
#' and the qlcMatrix package.
#' 
#' @param X A matrix
#' @param Y A matrix
#' @param cov return covariance matrix
#' @export
#' @author Michael Cysouw, Karsten Looschen
#' @importMethodsFrom Matrix colMeans
# covmat uses E[(X-muX)'(Y-muY)] = E[X'Y] - muX'muY
# with sample correction n/(n-1) this leads to cov = ( X'Y - n*muX'muY ) / (n-1)
#
# the sd in the case Y!=NULL uses E[X-mu]^2 = E[X^2]-mu^2
# with sample correction n/(n-1) this leads to sd^2 = ( X^2 - n*mu^2 ) / (n-1)
#
# Note that results larger than 1e4 x 1e4 will become very slow, because the resulting matrix is not sparse anymore.
corSparse <- function(X, Y = NULL, cov = FALSE) {
  X <- as(object = X, Class = "CsparseMatrix")
  n <- nrow(x = X)
  muX <- colMeans(x = X)

  if (!is.null(x = Y)) {
    if (nrow(x = X) != nrow(x = Y)) {
      stop("Matrices must contain the same number of rows")
    }
    Y <- as(object = Y, Class = "CsparseMatrix")
    muY <- colMeans(x = Y)
    covmat <- ( as.matrix(x = crossprod(x = X, y = Y)) - n * tcrossprod(x = muX, y = muY) ) / (n-1)
    sdvecX <- sqrt( (colSums(x = X^2) - n*muX^2) / (n-1) )
    sdvecY <- sqrt( (colSums(x = Y^2) - n*muY^2) / (n-1) )
    cormat <- covmat / tcrossprod(x = sdvecX, y = sdvecY)
  } else {		
    covmat <- ( as.matrix(crossprod(x = X)) - n * tcrossprod(x = muX) ) / (n-1)
    sdvec <- sqrt(x = diag(x = covmat))
    cormat <- covmat / tcrossprod(x = sdvec)
  }
  if (cov) {
    dimnames(x = covmat) <- NULL
    return(covmat)
  } else {
    dimnames(x = cormat) <- NULL
    return(cormat)
  }
}
