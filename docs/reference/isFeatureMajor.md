# Is optimized for feature access

Detect if matrix is designed for feature-wise access

## Usage

``` r
isFeatureMajor(x)
```

## Arguments

- x:

  file-backed `DelayedArray`

## Details

In R, check if the matrix has seed type `"CSR_H5ADMatrixSeed"`. This
corresponds to `"Compressed Sparse Column sparse matrix"` from
`anndata`. Note that a matrix optimized for accessing features (i.e.
genes) is coded as a CSC in python because the features are along
columns. But in R, this is refered to as CSR since features are along
rows. The underlying for matrix is the same, but calling `X[j,]` in R
returns all cells from gene `j`.

\# python code to get matrix type: \# import anndata \# ad =
anndata.read_h5ad(file) \# ad.X
