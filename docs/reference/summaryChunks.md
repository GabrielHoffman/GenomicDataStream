# Get information about chunks

Read through stream to get size of each chunk, IDs of variants and
samples

## Usage

``` r
summaryChunks(x)
```

## Arguments

- x:

  `GenomicDataStream`

## Examples

``` r
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

# initialize
obj <- GenomicDataStream(file, "DS", chunkSize = 3, initialize = TRUE)
summaryChunks(obj)
#> $intervals
#> [1] "1:10000-12000" "1:13000-15000" "1:16000-18000" "1:19000-19000"
#> 
#> $chunks
#> [1] 3 3 3 1
#> 
#> $sampleIDs
#>  [1] "I1"  "I2"  "I3"  "I4"  "I5"  "I6"  "I7"  "I8"  "I9"  "I10" "I11" "I12"
#> [13] "I13" "I14" "I15" "I16" "I17" "I18" "I19" "I20" "I21" "I22" "I23" "I24"
#> [25] "I25" "I26" "I27" "I28" "I29" "I30" "I31" "I32" "I33" "I34" "I35" "I36"
#> [37] "I37" "I38" "I39" "I40" "I41" "I42" "I43" "I44" "I45" "I46" "I47" "I48"
#> [49] "I49" "I50" "I51" "I52" "I53" "I54" "I55" "I56" "I57" "I58" "I59" "I60"
#> 
#> $variantIDs
#>  [1] "1:10000:C:A" "1:11000:T:C" "1:12000:T:C" "1:13000:T:C" "1:14000:G:C"
#>  [6] "1:15000:A:C" "1:16000:G:A" "1:17000:C:A" "1:18000:C:G" "1:19000:T:G"
#> 
#> $regions
#>  [1] "1:10000-10000" "1:11000-11000" "1:12000-12000" "1:13000-13000"
#>  [5] "1:14000-14000" "1:15000-15000" "1:16000-16000" "1:17000-17000"
#>  [9] "1:18000-18000" "1:19000-19000"
#> 
```
