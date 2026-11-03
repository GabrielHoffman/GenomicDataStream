

# write files
# FILE=~/prog/R-4.5.1/library/anndataR/extdata/example.h5ad
# OUT=~/workspace/repos/GenomicDataStream/inst/extdata

# recode_h5ad --input $FILE --output $OUT/example.h5ad

# recode_h5ad --input $FILE --output $OUT/example_csc.zarr --format CSC

# recode_h5ad --input $FILE --output $OUT/example_csr.zarr --format CSR


test_readZarr <- function(){

  library(SingleCellExperiment)
  library(GenomicDataStream)
  library(RUnit)

  # read zip'ed zarr to zarr
  getZarr <- function(zipfile, storeID){

    exdir <- file.path(tempdir(), "ZarrSparseMatrix-example")
    store <- file.path(exdir, storeID)
    if (!dir.exists(store)) {
      unzip(zipfile, exdir=exdir)
    }
    store
  }

  file <- system.file("extdata", "example.h5ad", package = "GenomicDataStream")
  sce <- readH5AD(file)

  file <- system.file("extdata", "example_csc.zarr.zip", package = "GenomicDataStream")
  loc2 = getZarr(file, "example_csc.zarr")
  sce2 = readZarr(loc2)

  file <- system.file("extdata", "example_csr.zarr.zip", package = "GenomicDataStream")
  loc3 = getZarr(file, "example_csr.zarr")
  sce3 = readZarr(loc3)

  counts(sce) = as.matrix(counts(sce))
  counts(sce2) = as.matrix(counts(sce2))
  counts(sce3) = as.matrix(counts(sce3))

  sce$libSize = colSums(counts(sce))
  names(sce$libSize) = c()

  metadata(sce) = list()
  metadata(sce2) = list()
  metadata(sce3) = list()

  checkEquals(colData(sce), colData(sce2), tol=1e-7)

  checkEquals(sce, sce2, tol=1e-7)
  checkEquals(sce, sce3, tol=1e-7)
}