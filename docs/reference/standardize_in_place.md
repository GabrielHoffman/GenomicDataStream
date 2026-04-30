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
#> [1,] -1.0489176 -0.9279937  0.8581271 -1.1534032
#> [2,]  0.9425843 -0.1310836  0.2400504  0.5293123
#> [3,]  0.1063333  1.0590773 -1.0981775  0.6240909
#> attr(,"scaled:center")
#> [1] 2.775558e-17 0.000000e+00 1.480297e-16 7.401487e-17
#> attr(,"scaled:scale")
#> [1] 1 1 1 1
```
