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
#> [1,]  0.9016198  0.3601163  0.6081774 -0.1704248
#> [2,] -1.0755589 -1.1301830 -1.1541418  1.0742607
#> [3,]  0.1739391  0.7700667  0.5459644 -0.9038359
#> attr(,"scaled:center")
#> [1] -2.775558e-17  0.000000e+00 -3.700743e-17 -3.700743e-17
#> attr(,"scaled:scale")
#> [1] 1 1 1 1
```
