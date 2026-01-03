# Plot PCAstream

Plot PCAstream

## Usage

``` r
# S4 method for class 'PCA,ANY'
plot(x, ..., main = "scree plot")
```

## Arguments

- x:

  `PCA` object

- ...:

  other arguments

- main:

  title

## Value

plot

## Examples

``` r
# Example: PCA on genotype data
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

obj <- GenomicDataStream(file, "DS", chunkSize = 3)

res <- PCAstream(obj, k=5, threads=1)
#> Read through...
#>  # features: 10
#>  # chunks: 2
#>  # PCs: 5
#>  # threads: 1
#> 
Epoch 0 / 7  
Epoch 1 / 7  
Epoch 2 / 7  
Epoch 3 / 7  
Epoch 4 / 7  
Epoch 5 / 7  
Epoch 6 / 7  
Epoch 7 / 7  
Final decompositions
Completed            

plot(res)
```
