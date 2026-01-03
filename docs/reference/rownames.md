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
#> NULL
```
