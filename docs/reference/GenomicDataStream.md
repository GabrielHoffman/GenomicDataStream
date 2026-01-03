# GenomicDataStream

Read genomic data files (VCF, BCF, BGEN, h5ad) into R/Rcpp in chunks for
analysis with Armadillo or Eigen libraries

Interface to GenomicDataStream C++ code

## Usage

``` r
GenomicDataStream(
  file,
  field = "",
  region = "",
  samples = "-",
  MAF = 0,
  minVariance = 0,
  chunkSize = 10000,
  missingToMean = TRUE,
  initialize = FALSE
)
```

## Arguments

- file:

  file in VCF/BCF/BGEN/PGEN format with index

- field:

  field of VCF/BCF to read. Ignored for other file types

- region:

  target in the format `chr2:1-12345`. Multiple regions can be separated
  by one of `",\n\t"`, for example `"chr2:1-12345, chr3:1000-8000"`.
  Setting region to `""` includes all variants

- samples:

  string of comma separated sample IDs to extract: `"ID1,ID2,ID3"`.
  `"-"` indicates all samples

- MAF:

  minor allele frequency filter applied to variants with max value
  \<= 2. Use `NaN` to retain all variants

- minVariance:

  minimum variance filter applied to variants with max value \> 2

- chunkSize:

  number of variants to return per chunk

- missingToMean:

  if true, set missing values to the mean dosage value. if false, set to
  `NaN`

- initialize:

  default `FALSE`. If `TRUE`, file info is read from path, otherwise
  store path until `GenomicDataStream` is initialized later

## Value

object of class `GenomicDataStream`

## Details

Variants are filtered using `MAF` if the max value is \<=2, or
`minVariance` otherwise

## Examples

``` r
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

# initialize
obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)

obj
#>       GenomicDataStream 
#> 
#>   file:          test.vcf.gz 
#>   initialized:   TRUE 
#>   stream type:   vcf.gz 
#>   field:         DS 
#>   region:         
#>   samples:       60 
#>   MAF:           0 
#>   minVar cutoff: 0 
#>   missingToMean: TRUE 
#>   chunkSize:     5 
#>   features read: 0 
#>   end of stream: FALSE 

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
```
