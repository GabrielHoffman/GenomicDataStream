# Set Chunk Size

Set chunk size for existing `GenomicDataStream`

## Usage

``` r
setChunkSize(x, chunkSize)
```

## Arguments

- x:

  `GenomicDataStream`

- chunkSize:

  positive integer

## Value

none

## Examples

``` r
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize=TRUE)

obj <- setChunkSize(obj, 200)
```
