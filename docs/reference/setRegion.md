# Set regions of GenomicDataStream

If `GenomicDataStream` is already initialized, set new query in C++
backend. Otherwise substitute `region` values

## Usage

``` r
setRegion(x, region)
```

## Arguments

- x:

  `GenomicDataStream`

- region:

  target in the format `chr2:1-12345`. Multiple regions can be separated
  by one of `",\n\t"`, for example `"chr2:1-12345, chr3:1000-8000"`.
  Setting region to `""` includes all variants

## Value

`GenomicDataStream` with region set

## Details

Set regions of GenomicDataStream

## Examples

``` r
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

obj <- GenomicDataStream(file, "DS", chunkSize = 5)

# by default, GenomicDataStream is not initialized
setRegion(obj, "1:10000-12000")
#>       GenomicDataStream 
#> 
#>   file:          test.vcf.gz 
#>   initialized:   FALSE 
```
