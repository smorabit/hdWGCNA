#' scale01
#'
#' Function to scale a numeric vector between 0 and 1.
#'
#' @param x numeric vector
#' @return A numeric vector the same length as x, rescaled to the range 0 to 1.
#' @keywords helper
#' @export
scale01 <- function(x) {
  y <- min(x)
  z <- max(x)
  (x - y) / (z - y)
}


#' umap_theme
#'
#' ggplot theme to remove axes etc.
#'
#' @return A ggplot2 theme object that removes axes, gridlines and backgrounds.
#' @keywords helper
#' @export
umap_theme <- function() {
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
}


#' shuffle_points
#'
#' Function to shuffle the rows of a dataframe.
#'
#' @param df dataframe
#' @return The input data frame with its rows in random order.
#' @keywords helper
#' @export
shuffle_points <- function(df) {
  df[sample(seq_len(nrow(df))), ]
}


#' Sparse matrix correlation
#'
#' Compute the Pearson correlation matrix between
#' columns of two sparse matrices.
#'
#' Originally from
#' \url{http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r}
#' and the qlcMatrix & Signac packages.
#'
#' @param X A matrix
#' @param Y A matrix
#' @param cov return covariance matrix
#' @return A dense matrix of Pearson correlations (or covariances if cov = TRUE) between the columns of X (and Y if provided), without dimnames.
#' @export
#' @author Michael Cysouw, Karsten Looschen
# covmat uses E[(X-muX)'(Y-muY)] = E[X'Y] - muX'muY
# with sample correction n/(n-1) this leads to cov = ( X'Y - n*muX'muY ) / (n-1)
#
# the sd in the case Y!=NULL uses E[X-mu]^2 = E[X^2]-mu^2
# with sample correction n/(n-1) this leads to sd^2 = ( X^2 - n*mu^2 ) / (n-1)
#
# Note that results larger than 1e4 x 1e4 will become very slow, because the resulting matrix is not sparse anymore.
corSparse <- function(X, Y = NULL, cov = FALSE) {
  X <- as(object = X, Class = 'CsparseMatrix')
  n <- nrow(x = X)
  muX <- colMeans(x = X)

  if (!is.null(x = Y)) {
    if (nrow(x = X) != nrow(x = Y)) {
      stop('Matrices must contain the same number of rows')
    }
    Y <- as(object = Y, Class = 'CsparseMatrix')
    muY <- colMeans(x = Y)
    covmat <- (as.matrix(x = crossprod(x = X, y = Y)) - n * tcrossprod(x = muX, y = muY)) / (n - 1)
    sdvecX <- sqrt((colSums(x = X^2) - n * muX^2) / (n - 1))
    sdvecY <- sqrt((colSums(x = Y^2) - n * muY^2) / (n - 1))
    cormat <- covmat / tcrossprod(x = sdvecX, y = sdvecY)
  } else {
    covmat <- (as.matrix(crossprod(x = X)) - n * tcrossprod(x = muX)) / (n - 1)
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


#' CheckSeurat5
#'
#' Function to check if the installed Seurat is version 5 or newer
#'
#' @return TRUE if Seurat >= 5.0.0 is installed
#' @export
CheckSeurat5 <- function() {
  utils::packageVersion('Seurat') >= '5.0.0'
}


# Bioconductor packages among hdWGCNA's suggested dependencies (for install hints)
bioc_suggests <- c(
  'BiocGenerics', 'DESeq2', 'ensembldb', 'GeneOverlap', 'GenomeInfoDb', 'GenomicRanges',
  'IRanges', 'motifmatchr', 'SingleCellExperiment', 'UCell'
)


# GitHub-only packages used by optional features (not installable from CRAN/Bioconductor)
github_suggests <- c(ggradar = 'ricardo-bion/ggradar')


# Stop with an informative message if packages needed by an optional feature are missing.
# Use at the top of any function that calls a package listed in Suggests.
check_suggests <- function(pkgs, feature = NULL) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    return(invisible(TRUE))
  }

  hints <- character(0)
  cran <- setdiff(missing, c(bioc_suggests, names(github_suggests)))
  bioc <- intersect(missing, bioc_suggests)
  gh <- intersect(missing, names(github_suggests))
  if (length(cran) > 0) {
    hints <- c(hints, paste0("install.packages(c('", paste(cran, collapse = "', '"), "'))"))
  }
  if (length(bioc) > 0) {
    hints <- c(hints, paste0("BiocManager::install(c('", paste(bioc, collapse = "', '"), "'))"))
  }
  if (length(gh) > 0) {
    hints <- c(hints, paste0("remotes::install_github('", github_suggests[gh], "')"))
  }

  stop(
    paste0(
      ifelse(length(missing) > 1, 'Packages ', 'Package '),
      paste0("'", missing, "'", collapse = ', '),
      ifelse(length(missing) > 1, ' are', ' is'),
      ' required',
      ifelse(is.null(feature), '', paste0(' for ', feature)),
      '. Install with:\n  ',
      paste(hints, collapse = '\n  ')
    ),
    call. = FALSE
  )
}


# TRUE if x is a base matrix with numeric values (replaces tester::is_numeric_matrix)
is_numeric_matrix <- function(x) {
  is.matrix(x) && is.numeric(x)
}


# TRUE if x is a data.frame whose columns are all numeric
is_numeric_dataframe <- function(x) {
  is.data.frame(x) && all(vapply(x, is.numeric, logical(1)))
}


# Stop with an informative message if a dimensionality reduction is missing.
check_reduction <- function(seurat_obj, reduction) {
  if (!(reduction %in% names(seurat_obj@reductions))) {
    stop(paste0(
      'Reduction ', reduction, ' not found in seurat_obj. Available reductions: ',
      paste0(names(seurat_obj@reductions), collapse = ', '), '.'
    ), call. = FALSE)
  }
  invisible(TRUE)
}


# Alternating horizontal background stripes for plots with a discrete y axis
# (e.g. one row per module). Same approach as ggforestplot::geom_stripes (MIT licence,
# NightingaleHealth), reimplemented here so hdWGCNA doesn't depend on a GitHub-only package.
geom_stripes <- function(
  mapping = NULL,
  data = NULL,
  stat = 'identity',
  position = 'identity',
  ...,
  show.legend = NA,
  inherit.aes = TRUE
) {
  ggplot2::layer(
    geom = GeomStripes,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(...)
  )
}


GeomStripes <- ggplot2::ggproto(
  'GeomStripes',
  ggplot2::Geom,
  required_aes = c('y'),
  default_aes = ggplot2::aes(
    xmin = -Inf, xmax = Inf,
    odd = '#22222222', even = '#00000000',
    alpha = NA, colour = NA, linetype = 'solid', linewidth = 0.5
  ),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(data, panel_params, coord) {
    y <- round(data$y)
    rects <- unique(data.frame(
      xmin = data$xmin, xmax = data$xmax,
      ymin = y - 0.5, ymax = y + 0.5,
      odd = data$odd, even = data$even,
      alpha = data$alpha, colour = data$colour,
      linetype = data$linetype, linewidth = data$linewidth,
      PANEL = data$PANEL, group = -1
    ))
    rects <- rects[order(rects$ymin), ]
    rects$fill <- ifelse(seq_len(nrow(rects)) %% 2 == 1, rects$odd, rects$even)
    rects$odd <- NULL
    rects$even <- NULL
    ggplot2::GeomRect$draw_panel(rects, panel_params, coord)
  }
)


# Recover the absolute TOM path from paths stored by hdWGCNA <= 0.4.13 as
# paste0(getwd(), '/', <absolute path>). Returns NA if the path doesn't look like that.
repair_tom_path <- function(path) {
  # Windows: the last drive letter starts the real path, e.g. 'C:/x/D:/data/TOM.rda'
  drives <- gregexpr('[A-Za-z]:[/\\\\]', path)[[1]]
  if (length(drives) > 1 && drives[1] != -1) {
    return(substring(path, drives[length(drives)]))
  }
  # Unix: '/home/me//data/TOM.rda'
  if (grepl('//', path, fixed = TRUE)) {
    return(sub('^.*//', '/', path))
  }
  NA_character_
}
