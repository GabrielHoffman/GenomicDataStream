# Reinitialize GenomicDataStream

Read file info from path to initialise stream

## Usage

``` r
reinitializeStream(x, region = NULL)
```

## Arguments

- x:

  `GenomicDataStream`

- region:

  new set of region, as a string

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

reinitializeStream(obj)
#>       GenomicDataStream 
#> 
#>   file:          test.vcf.gz 
#>   initialized:   TRUE 
#>   stream type:   vcf.gz 
#>   field:         DS 
#>   region:         
#>   samples:       60 
#>   minVar cutoff: 0 
#>   missingToMean: TRUE 
#>   chunkSize:     5 
#>   features read: 0 
#>   end of stream: FALSE 
```
