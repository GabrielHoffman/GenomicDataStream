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
