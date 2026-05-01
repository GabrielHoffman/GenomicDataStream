# Evaluate performance of PC estimates

Evaluate performance of PC estimates compared to true PC values

## Usage

``` r
perfMetric(U, U_est, k = ncol(U_est), metric = c("MEV", "minSSE"))
```

## Arguments

- U:

  true eigen values

- U_est:

  estimated eigen-values

- k:

  number of PCs to evaluate

- metric:

  evaluate the accuracy of the estimated PCs compared to the true PCs
  using mean explained variance (`"MEV"`) or minimum of sum of squared
  errors (`"minSSE"`)

## Value

performance metric

## Details

See performance metrics described by Li, et al. (2023)

## References

- Li, Z., Meisner, J., & Albrechtsen, A. (2023). Fast and accurate
  out-of-core PCA framework for large scale biobank data. Genome
  Research, 33(9), 1599-1608.
  [doi:10.1101/gr.277525.122](https://doi.org/10.1101/gr.277525.122) .

## Examples

``` r
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
X <- hilbert(9)[, 1:6]
k <- 4

# Compute SVD using two methods
dcmp1 <- svd(X, k, k)
dcmp2 <- dashSVD( X, k=k)

# Mean variance explained is 1
perfMetric(dcmp1$u, dcmp2$u, metric = "MEV")
#> [1] 1

# minimum of sum of squared errors zero
perfMetric(dcmp1$u, dcmp2$u, metric = "minSSE")
#> [1] 1.100796e-12
```
