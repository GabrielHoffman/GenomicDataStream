# Standardize matrix columns in place

Standardize mean and variance of matrix columns in place

## Usage

``` r
standardize_in_place(X, center = TRUE, scale = TRUE)
```

## Arguments

- X:

  matrix

- center:

  boolean, TRUE indices center columns

- scale:

  boolean, TRUE indices scale columns

## Value

none, matrix is standardized in place

## Examples

``` r
X <- matrix(runif(12), 3, 4)
standardize_in_place(X)

scale(X)
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  0.4101263 -0.3965980 -0.4988723  0.9336013
#> [2,] -1.1398612  1.1374649  1.1512925 -1.0552660
#> [3,]  0.7297349 -0.7408669 -0.6524202  0.1216647
#> attr(,"scaled:center")
#> [1] 3.700743e-17 1.110223e-16 1.110223e-16 2.081668e-16
#> attr(,"scaled:scale")
#> [1] 1 1 1 1
```
