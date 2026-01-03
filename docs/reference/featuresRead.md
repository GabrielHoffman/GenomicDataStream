# Get number of features read from GenomicDataStream

Get number of total features read from GenomicDataStream

## Usage

``` r
featuresRead(x)
```

## Arguments

- x:

  `GenomicDataStream`

## Value

total number of features read from the stream

## Examples

``` r
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

# initialize
obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)

# loop until break
while (1) {
  # get data chunk
  # data$X matrix with features as columns
  # data$info information about each feature as rows
  dat <- getNextChunk(obj)

  if (atEndOfStream(obj)) break

  print(dat$info)
}
#>   CHROM   POS          ID A1 A2
#> 1     1 10000 1:10000:C:A  C  A
#> 2     1 11000 1:11000:T:C  T  C
#> 3     1 12000 1:12000:T:C  T  C
#> 4     1 13000 1:13000:T:C  T  C
#> 5     1 14000 1:14000:G:C  G  C
#>   CHROM   POS          ID A1 A2
#> 1     1 15000 1:15000:A:C  A  C
#> 2     1 16000 1:16000:G:A  G  A
#> 3     1 17000 1:17000:C:A  C  A
#> 4     1 18000 1:18000:C:G  C  G
#> 5     1 19000 1:19000:T:G  T  G

featuresRead(obj)
#> [1] 10
```
