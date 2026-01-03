# Normalize principal components

Modify signs of principal components so diagonal are always positive

## Usage

``` r
normPC(U)
```

## Arguments

- U:

  matrix of principal components

## Value

matrix of principal components with positive diagonal values

## Examples

``` r
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
 X <- hilbert(9)[, 1:6]
k <- 4

dcmp <- svd( scale(X), k, k)

normPC( dcmp$u ) 
#>               [,1]         [,2]        [,3]        [,4]
#>  [1,]  0.733029934 -0.560955931  0.18659115  0.04466128
#>  [2,]  0.346439868  0.398217022 -0.65845777 -0.39266193
#>  [3,]  0.133584689  0.436130386  0.02675310  0.57827877
#>  [4,] -0.008363917  0.311719143  0.34136688  0.23654274
#>  [5,] -0.111482887  0.157448097  0.37596273 -0.19400608
#>  [6,] -0.190368173  0.006937458  0.25683645 -0.38066482
#>  [7,] -0.252905441 -0.130842047  0.06005063 -0.30216922
#>  [8,] -0.303813760 -0.254331745 -0.17222059 -0.01357099
#>  [9,] -0.346120314 -0.364322384 -0.41688258  0.42359025
```
