# Get genomic intervals

Get genomic intervals chopping up regions into smaller chunks

## Usage

``` r
chopChroms(x, nchunks = 10)
```

## Arguments

- x:

  `data.frame` with columns chrom, start, end

- nchunks:

  number of chunks to create

## Value

vector genome intervals

## Examples

``` r
file <- system.file("extdata", "test.bed", package = "GenomicDataStream")

# initialize
obj <- GenomicDataStream(file, chunkSize = 5, initialize = TRUE)
df <- getChromRanges(obj)
chopChroms( df )
#>  [1] "1:10000-10900" "1:10901-11800" "1:11801-12700" "1:12701-13600"
#>  [5] "1:13601-14500" "1:14501-15400" "1:15401-16300" "1:16301-17200"
#>  [9] "1:17201-18100" "1:18101-19000"
```
