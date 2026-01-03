# Initialize GenomicDataStream

Read file info from path to initialise stream. If already initialized,
return to the beginning of the stream

## Usage

``` r
initializeStream(x, region = NULL)
```

## Arguments

- x:

  `GenomicDataStream`

- region:

  target in the format `chr2:1-12345`. Multiple regions can be separated
  by one of `",\n\t"`, for example `"chr2:1-12345, chr3:1000-8000"`.
  Setting region to `""` includes all variants

## Value

initialized `GenomicDataStream`

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
