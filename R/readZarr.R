


#' Read Zarr as SingleCellExperiment
#'
#' Read Zarr as SingleCellExperiment where counts is a file-backed DelayedArray
#' 
#' @param file Zarr file 
#' @param target (default: \code{"X"}) name of the matrix to read 
#' @param ondisk if \code{TRUE} (default), only stream count data into memory when needed.  If \code{FALSE}, read count data into memory now as a \code{sparseMatrix}
#' @param verbose print messages
#' 
#' @return \code{SingleCellExperiment}
#' 
#' @examples
#' # read zip'ed zarr to zarr
#' getZarr <- function(zipfile, storeID){
#' 
#'   exdir <- file.path(tempdir(), "ZarrSparseMatrix-example")
#'   store <- file.path(exdir, storeID)
#'   if (!dir.exists(store)) {
#'     unzip(zipfile, exdir=exdir)
#'   }
#'   store
#' }
#' 
#' loc <- getZarr(system.file("extdata", "example_csc.zarr.zip", package = "GenomicDataStream"), "example_csc.zarr")
#' 
#' sce <- readZarr(loc)
#' 
#' sce
# 
#' @importFrom ZarrSparseMatrix ZarrSparseMatrix
#' @importFrom anndataR read_zarr
#' @export
readZarr <- function(file, target="X", ondisk = TRUE, verbose=FALSE){

  if( ! file.exists(file) ){
    txt <- paste("File does not exist:", file)
    stop(txt)
  }

  if( verbose ){
    message("Reading counts...")
  }

  # Read data as Delayed/HDF5-backed matrix
  # if layer = NULL, read X.  Otherwise use layer name
  # returns observations matrix (genes x cells)
  counts <- tryCatch({
    t(ZarrSparseMatrix(file, target))
    }, 
    error = function(e){
      e$message <- paste0("Error reading file: ", conditionMessage(e))
      e
      })

  if( verbose ){
    axis <- ifelse( isFeatureMajor(counts), "genes", "cells")
    txt <- paste("  File optimized for accessing:", axis)
    message(txt)
  }

  if( ! ondisk ){
    if( verbose ){
      message("  Converting to in-memory sparse matrix...")
    }
    counts <- as(counts, "dgCMatrix")
  }

  if( verbose ){
    message("Reading supporting data...")
  }

  # Read Zarr, file-backed until data access 
  ad <- read_zarr(file, "ZarrAnnData")

  if( verbose ){
    message("  Reading colData (/obs)...")
  }
  colData <- data.frame(ad$obs)

  if( verbose ){
    message("  Reading rowData (/var)...")
  }
  rowData <- data.frame(ad$var)

  dimnames(counts) <- list(rownames(rowData), rownames(colData))

  # Wrap as a SingleCellExperiment
  # only X is file-backed
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = colData,
    rowData = rowData,
    metadata = ad$uns)

  sce
}

