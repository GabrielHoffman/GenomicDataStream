# Get status of GenomicDataStream

If `initialized`, return `TRUE`, else `FALSE`

## Usage

``` r
isInitialized(x)
```

## Arguments

- x:

  `GenomicDataStream`

## Value

initialization status

## Examples

``` r
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

obj <- GenomicDataStream(file, "DS", chunkSize = 5)

# by default, GenomicDataStream is not initialized
isInitialized(obj)
#> [1] FALSE

# initialize
obj <- initializeStream(obj)

isInitialized(obj)
#> [1] TRUE
```
