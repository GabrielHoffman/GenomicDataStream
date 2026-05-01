# Read H5AD as SingleCellExperiment

Read H5AD as SingleCellExperiment where counts is a file-backed
DelayedArray

## Usage

``` r
readH5AD(file, layer = NULL, ondisk = TRUE, verbose = FALSE, raw = FALSE)
```

## Arguments

- file:

  H5AD file

- layer:

  `NULL` (the default) or the name of a matrix in the `/layers` group.
  By default (i.e. when `layer` is not specified) returns the central
  matrix (`X`).

- ondisk:

  if `TRUE` (default), only stream count data into memory when needed.
  If `FALSE`, read count data into memory now as a `sparseMatrix`

- verbose:

  print messages

- raw:

  if `TRUE`, read counts from `/raw/X`. Cannot be used with `layer`.

## Value

`SingleCellExperiment`

## Details

Uses
[`HDF5Array::H5ADMatrix()`](https://rdrr.io/pkg/HDF5Array/man/H5ADMatrix-class.html)
to read counts as a file-backed DelayedArray, and
[`anndataR::read_h5ad()`](https://anndataR.scverse.org/reference/read_h5ad.html)
to read all other data from H5AD.

## Examples

``` r
file <- system.file("extdata", "example.h5ad", package = "anndataR")
sce <- readH5AD(file)
```
