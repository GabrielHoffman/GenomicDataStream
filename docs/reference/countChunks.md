# Get chunk sizes

Get array of chunk sizes

## Usage

``` r
countChunks(x)
```

## Arguments

- x:

  `GenomicDataStream`

## Value

array of chunk sizes

## Examples

``` r
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

# initialize
obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)

countChunks(obj)
#> [1] 5 5
```
