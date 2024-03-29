% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{corSparse}
\alias{corSparse}
\title{Sparse matrix correlation}
\usage{
corSparse(X, Y = NULL, cov = FALSE)
}
\arguments{
\item{X}{A matrix}

\item{Y}{A matrix}

\item{cov}{return covariance matrix}
}
\description{
Compute the Pearson correlation matrix between
columns of two sparse matrices.
}
\details{
Originally from
\url{http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r}
and the qlcMatrix & Signac packages.
}
\author{
Michael Cysouw, Karsten Looschen
}
