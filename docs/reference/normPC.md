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

dcmp <- svd(X, k, k)

normPC( dcmp$u ) 
#>            [,1]       [,2]        [,3]        [,4]
#>  [1,] 0.7244999 -0.6265620 -0.27350003  0.08526902
#>  [2,] 0.4281556  0.1298781  0.64293597 -0.55047428
#>  [3,] 0.3121985  0.2803679  0.33633240  0.31418014
#>  [4,] 0.2478932  0.3141885  0.06931246  0.44667149
#>  [5,] 0.2063780  0.3140734 -0.10786005  0.30241655
#>  [6,] 0.1771408  0.3026808 -0.22105904  0.09041508
#>  [7,] 0.1553452  0.2877310 -0.29280775 -0.11551327
#>  [8,] 0.1384280  0.2721599 -0.33783778 -0.29312535
#>  [9,] 0.1248940  0.2571250 -0.36542543 -0.43884649
```
