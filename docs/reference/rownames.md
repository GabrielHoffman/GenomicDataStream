# Get rownames

BGEN uses sample order from the query, but VCF/BCF/PGEN uses order in
file

## Usage

``` r
# S4 method for class 'GenomicDataStream'
rownames(x, do.NULL = TRUE, prefix = "row")
```

## Arguments

- x:

  `GenomicDataStream`

- do.NULL:

  not used

- prefix:

  not used

## Value

array of string names

## Details

Get rownames (i.e. sample names) in order that the samples are extracted

## Examples

``` r
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

obj <- GenomicDataStream(file, "DS", init=TRUE)

rownames(obj)
#>  [1] "I1"  "I2"  "I3"  "I4"  "I5"  "I6"  "I7"  "I8"  "I9"  "I10" "I11" "I12"
#> [13] "I13" "I14" "I15" "I16" "I17" "I18" "I19" "I20" "I21" "I22" "I23" "I24"
#> [25] "I25" "I26" "I27" "I28" "I29" "I30" "I31" "I32" "I33" "I34" "I35" "I36"
#> [37] "I37" "I38" "I39" "I40" "I41" "I42" "I43" "I44" "I45" "I46" "I47" "I48"
#> [49] "I49" "I50" "I51" "I52" "I53" "I54" "I55" "I56" "I57" "I58" "I59" "I60"
```
