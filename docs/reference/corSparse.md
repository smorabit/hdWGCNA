# Sparse matrix correlation

Compute the Pearson correlation matrix between columns of two sparse
matrices.

## Usage

``` r
corSparse(X, Y = NULL, cov = FALSE)
```

## Arguments

- X:

  A matrix

- Y:

  A matrix

- cov:

  return covariance matrix

## Details

Originally from
<http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r>
and the qlcMatrix & Signac packages.

## Author

Michael Cysouw, Karsten Looschen
