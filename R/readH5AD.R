
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
  is(seed(x), "CSR_H5ADMatrixSeed") || is(seed(x), "CSR_H5SparseMatrixSeed")
}


#' Read H5AD as SingleCellExperiment
#'
#' Read H5AD as SingleCellExperiment where counts is a file-backed DelayedArray
#' 
#' @param file H5AD file 
#' @param layer \code{NULL} (the default) or the name of a matrix in the \code{/layers} group. By default (i.e. when \code{layer} is not specified) returns the central matrix (\code{X}).
#' @param ondisk if \code{TRUE} (default), only stream count data into memory when needed.  If \code{FALSE}, read count data into memory now as a \code{sparseMatrix}
#' @param verbose print messages
#' @param raw if \code{TRUE}, read counts from \code{/raw/X}. Cannot be used with \code{layer}.
#' 
#' @return \code{SingleCellExperiment}
#' 
#' @examples
#' file <- system.file("extdata", "example.h5ad", package = "anndataR")
#' sce <- readH5AD(file)
# 
#' @details Uses \code{HDF5Array::H5ADMatrix()} to read counts as a file-backed DelayedArray, and \code{anndataR::read_h5ad()} to read all other data from H5AD.
#' 
#' @importFrom HDF5Array H5ADMatrix
#' @importFrom HDF5Array H5SparseMatrix
#' @importFrom anndataR read_h5ad
#' @importFrom methods as
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom rhdf5filters hdf5_plugin_path
#' @importFrom withr with_envvar
#' @export
readH5AD <- function(file, layer=NULL, ondisk = TRUE, verbose=FALSE, raw=FALSE){

  if( ! file.exists(file) ){
    txt <- paste("File does not exist:", file)
    stop(txt)
  }

  if( ! isTRUE(raw) && ! identical(raw, FALSE) ){
    stop("'raw' must be TRUE or FALSE")
  }

  if( isTRUE(raw) && ! is.null(layer) ){
    stop("'raw=TRUE' cannot be used with 'layer'")
  }
  
  # Load hdf5 compression plugins: not allowed in BioC
  # Sys.setenv("HDF5_PLUGIN_PATH" = hdf5_plugin_path())

  if( verbose ){
    message("Reading counts...")
  }

  if( isTRUE(raw)  ){
    # check that raw/X group exists
    with_envvar(
        new = c(HDF5_PLUGIN_PATH = hdf5_plugin_path()),
        code = .check_group(file, "raw/X"))
  }

  # Read data as Delayed/HDF5-backed matrix
  # if layer = NULL, read X.  Otherwise use layer name
  # returns observations matrix (genes x cells)
  tryCatch({

    if( isTRUE(raw) ){
      # set HDF5_PLUGIN_PATH only within the local env
      # without using Sys.setenv
      counts <- with_envvar(
        new = c(HDF5_PLUGIN_PATH = hdf5_plugin_path()),
        code = H5SparseMatrix(file, "raw/X"))
    }else{
      counts <- with_envvar(
        new = c(HDF5_PLUGIN_PATH = hdf5_plugin_path()),
        code = H5ADMatrix(file, layer=layer)) 
    }
    }, 
    error = function(e){
      e$message <- paste0("Error reading file: ", conditionMessage(e), ". Issue with compression plugin?")
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

  # Read H5AD, file-backed until data access 
  ad <- read_h5ad(file, "HDF5AnnData")

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


# Adapted from HDF5Array:::.check_group()
#' @importFrom h5mread h5exists h5isdataset h5isgroup
#' @importFrom S4Vectors wmsg 
.check_group <- function (filepath, group){
    if (!h5exists(filepath, group)) 
        stop(wmsg("HDF5 group \"", group, "\" does not exist ", 
            "in this HDF5 file"))
    if (h5isdataset(filepath, group)) {
        is_h5ad_X_or_layer <- group == "/X" || substr(group, 
            1L, 8L) == "/layers/"
        msg1 <- c("\"", group, "\" is an HDF5 dataset, not an HDF5 group, ", 
            "so it looks like the matrix that you are trying to ", 
            "access is not stored in a sparse format. Please ", 
            "consider using the ")
        if (is_h5ad_X_or_layer) {
            msg2 <- c("H5ADMatrix() constructor if you are trying ", 
                "to access the central matrix of an h5ad file. ", 
                "Otherwise, use the HDF5Array() constructor.")
        }
        else {
            msg2 <- "HDF5Array() constructor to access this dataset."
        }
        stop(wmsg(msg1, msg2))
    }
    if (!h5isgroup(filepath, group)) 
        stop(wmsg("HDF5 object \"", group, "\" is not a group"))
}
