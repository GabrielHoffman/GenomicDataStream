# Collapse array of regions

Collapse array of regions into one interval per chromosome

## Usage

``` r
collapseRegions(regions)
```

## Arguments

- regions:

  array of interval strings

## Value

`data.frame` one interval per chromosome spanning all given regions

## Examples

``` r
regions <- c( "chr2:3-5", "chr1:4-5", "chr1:1-3")
collapseRegions( regions )
#>   chrom start end
#> 1  chr1     1   5
#> 2  chr2     3   5
```
