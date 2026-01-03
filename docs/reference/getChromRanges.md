# Get ranges for each chromosome

Get max and max postion for each chromosome

## Usage

``` r
getChromRanges(x)
```

## Arguments

- x:

  `GenomicDataStream`

## Value

`data.frame` of chrom, start, end

## Examples

``` r
file <- system.file("extdata", "test.bed", package = "GenomicDataStream")

# initialize
obj <- GenomicDataStream(file, chunkSize = 5, initialize = TRUE)
getChromRanges( obj )
#>   chrom start   end
#> 1     1 10000 19000
```
