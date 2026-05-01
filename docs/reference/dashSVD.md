# Dynamic Shifts Based Randomized SVD

DashSVD for fast randomized SVD

## Usage

``` r
dashSVD(x, k, p = 3, s = 20, rand = 1, byrow = FALSE, verbose = FALSE)

# S4 method for class 'ANY'
dashSVD(x, k, p = 3, s = 20, rand = 1, byrow = FALSE, verbose = FALSE)
```

## Arguments

- x:

  data matrix

- k:

  specifies the target rank of the low-rank decomposition. \\k\\ should
  satisfy \\k \<\< min(m,n)\\.

- p:

  number of additional power iterations (by default \\p=7\\).

- s:

  oversampling parameter (by default \\s=20\\).

- rand:

  randomization method. 1 indicates Gaussian sampling and 2 indicates
  uniform sampling

- byrow:

  default `FALSE`. Set to true if accessing rows is faster that
  accessing columns

- verbose:

  print messages

## Value

`PCA` class storing `d`, `u` and `v`

## References

Feng, X., Yu, W., Xie, Y. and Tang, J., 2024. Algorithm 1043: Faster
randomized svd with dynamic shifts. ACM Transactions on Mathematical
Software, 50(2), pp.1-27.

## Author

Zilong Li <zilong.dk@gmail.com>, Gabriel Hoffman

## Examples

``` r

hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
X <- hilbert(9)[, 1:6]

dcmp1 <- svd(X)
dcmp2 <- dashSVD(X, k=5)

# compare singular values
dcmp1$d
#> [1] 1.668433e+00 2.773727e-01 2.223722e-02 1.084693e-03 3.243788e-05
#> [6] 5.234864e-07
dcmp2$d
#> [1] 1.668433194 0.277372680 0.022237220 0.001084693 0.000000000
```
