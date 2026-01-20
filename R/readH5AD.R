
#' Is optimized for feature access
#'
#' Detect if matrix is designed for feature-wise access
#'
#' @param x file-backed \code{DelayedArray}
#' 
#' @details In R, check if the matrix has seed type \code{"CSR_H5ADMatrixSeed"}.  This corresponds to \code{"Compressed Sparse Column sparse matrix"} from \code{anndata}.  Note that a matrix optimized for accessing features (i.e. genes) is coded as a CSC in python because the features are along columns.  But in R, this is refered to as CSR since features are along rows.  The underlying for matrix is the same, but calling \code{X[j,]} in R returns all cells from gene \code{j}.
#' 
#' # python code to get matrix type:
#' # import anndata
#' # ad = anndata.read_h5ad(file)
#' # ad.X
#' @importFrom DelayedArray seed
#' @keywords internal
isFeatureMajor <- function(x){
  is(seed(x), "CSR_H5ADMatrixSeed")
}


#' Read H5AD as SingleCellExperiment
#'
#' Read H5AD as SingleCellExperiment where counts is a file-backed DelayedArray
#' 
#' @param file H5AD file 
#' @param layer \code{NULL} (the default) or the name of a matrix in the \code{/layers} group. By default (i.e. when \code{layer} is not specified) returns the central matrix (\code{X}).
#' @param verbose print messages
#' 
#' @details Uses \code{HDF5Array::H5ADMatrix()} to read counts as a file-backed DelayedArray, and \code{anndataR::read_h5ad()} to read all other data from H5AD.
#' 
#' @importFrom HDF5Array H5ADMatrix
#' @importFrom anndataR read_h5ad
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom rhdf5filters hdf5_plugin_path
#' @export
readH5AD <- function(file, layer=NULL, verbose=FALSE){
  
  # Load hdf5 compression plugins
  Sys.setenv("HDF5_PLUGIN_PATH" = hdf5_plugin_path())

  # Read data as Delayed/HDF5-backed matrix
  # if layer = NULL, read X.  Otherwise use layer name
  # returns observations matrix (genes x cells)
  tryCatch({
    counts <- H5ADMatrix(file, layer=layer) 
    }, 
    error = function(e){
      stop("Error reading file. Issue with compression plugin?")
      })

  if( verbose ){
    axis <- ifelse( isFeatureMajor(counts), "genes", "cells")
    txt <- paste("File optimized for accessing:", axis)
    message(txt)
  }

  # Read H5AD, file-backed until data access 
  ad <- read_h5ad(file, "HDF5AnnData")

  # Wrap as a SingleCellExperiment
  # only X is file-backed
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = data.frame(ad$obs),
    rowData = data.frame(ad$var),
    metadata = ad$uns)

  sce
}
