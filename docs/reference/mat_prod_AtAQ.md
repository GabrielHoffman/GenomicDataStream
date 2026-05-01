# Evaluate `crossprod(A, A %*% Q)`

When `A` is a `matrix`, use standard evaluation. When `A` is
`ResidualMatrixGLM`, use block-wise evaluation

## Usage

``` r
mat_prod_AtAQ(A, Q)

# S4 method for class 'ANY'
mat_prod_AtAQ(A, Q)
```

## Arguments

- A:

  any valid matrix type or `ResidualMatrixGLM`

- Q:

  any valid matrix type

## Value

`crossprod(A, A %*% Q)`

## Examples

``` r
A <- matrix(rnorm(4), 2, 2)
Q <- matrix(rnorm(2), 2, 1)

crossprod(A, A %*% Q)
#>          [,1]
#> [1,] 2.735281
#> [2,] 1.929868
mat_prod_AtAQ(A, Q)
#>          [,1]
#> [1,] 2.735281
#> [2,] 1.929868
```
