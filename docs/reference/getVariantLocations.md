# Get location of each feature

Get location of each feature in GenomicDataStream

## Usage

``` r
getVariantLocations(x)
```

## Arguments

- x:

  `GenomicDataStream`

## Value

get data chunk as `list` with entries `X` storing a matrix with features
as columns, and `info` storing information about each feature as rows

## Examples

``` r
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

# initialize
obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)

getVariantLocations(obj)
#>  [1] "1:10000-10000" "1:11000-11000" "1:12000-12000" "1:13000-13000"
#>  [5] "1:14000-14000" "1:15000-15000" "1:16000-16000" "1:17000-17000"
#>  [9] "1:18000-18000" "1:19000-19000"
```
