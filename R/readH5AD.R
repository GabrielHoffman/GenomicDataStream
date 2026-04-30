
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
#' @details Uses \code{HDF5Array::H5ADMatrix()} to read counts as a file-backed DelayedArray, and \code{anndataR::read_h5ad()} to read all other data from H5AD.
#' 
#' @importFrom HDF5Array H5ADMatrix
#' @importFrom HDF5Array H5SparseMatrix
#' @importFrom anndataR read_h5ad
#' @importFrom methods as
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom rhdf5filters hdf5_plugin_path
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
  
  # Load hdf5 compression plugins
  Sys.setenv("HDF5_PLUGIN_PATH" = hdf5_plugin_path())

  if( verbose ){
    message("Reading counts...")
  }

  if( isTRUE(raw) && ! .h5ad_path_exists(file, "raw/X") ){
    stop("Cannot read raw counts: HDF5 object /raw/X does not exist in this H5AD file")
  }

  # Read data as Delayed/HDF5-backed matrix
  # if layer = NULL, read X.  Otherwise use layer name
  # returns observations matrix (genes x cells)
  tryCatch({
    if( isTRUE(raw) ){
      counts <- H5SparseMatrix(file, "raw/X")
    }else{
      counts <- H5ADMatrix(file, layer=layer) 
    }
    }, 
    error = function(e){
      stop("Error reading file: ", conditionMessage(e), ". Issue with compression plugin?")
      })

  if( verbose ){
    axis <- ifelse( isFeatureMajor(counts), "genes", "cells")
    txt <- paste("File optimized for accessing:", axis)
    message(txt)
  }

  if( ! ondisk ){
    counts <- as(counts, "dgCMatrix")
  }

  if( verbose ){
    message("Reading colData, etc...")
  }

  # Read H5AD, file-backed until data access 
  ad <- read_h5ad(file, "HDF5AnnData")

  colData <- data.frame(ad$obs)
  if( isTRUE(raw) && .h5ad_path_exists(file, "raw/var") ){
    rowData <- .read_h5ad_dataframe(file, "raw/var")
  }else{
    rowData <- data.frame(ad$var)
  }

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


#' @importFrom rhdf5 h5read h5readAttributes h5ls
.read_h5ad_vector <- function(file, path){
  otype <- .h5ad_object_type(file, path)

  if( identical(otype, "H5I_DATASET") ){
    return(suppressWarnings(h5read(file, path)))
  }

  if( identical(otype, "H5I_GROUP") ){
    attrs <- suppressWarnings(h5readAttributes(file, path))
    encoding_type <- attrs[["encoding-type"]]
    if( identical(as.character(encoding_type), "categorical") ){
      codes <- suppressWarnings(h5read(file, file.path(path, "codes")))
      categories <- suppressWarnings(h5read(file, file.path(path, "categories")))
      ans <- rep(NA_character_, length(codes))
      keep <- !is.na(codes) & codes >= 0
      ans[keep] <- as.character(categories[codes[keep] + 1L])
      return(ans)
    }
  }

  stop(sprintf("Unsupported H5AD dataframe column: %s", path))
}

.h5ad_object_type <- function(file, path){
  path <- sub("^/", "", path)
  listing <- suppressWarnings(h5ls(file, recursive = TRUE))
  object_paths <- file.path(sub("^/", "", listing$group), listing$name)
  object_paths <- sub("^\\./", "", object_paths)
  object_paths <- sub("^/", "", object_paths)
  hit <- which(object_paths == path)
  if( ! length(hit) ){
    return(NA_character_)
  }
  as.character(listing$otype[hit[[1]]])
}

.h5ad_path_exists <- function(file, path){
  ! is.na(.h5ad_object_type(file, path))
}

.read_h5ad_dataframe <- function(file, group){
  attrs <- suppressWarnings(h5readAttributes(file, group))
  index_col <- attrs[["_index"]]
  columns <- attrs[["column-order"]]
  fields <- unique(c(index_col, columns))
  fields <- fields[nzchar(fields)]

  values <- lapply(fields, function(field){
    .read_h5ad_vector(file, file.path(group, field))
  })
  names(values) <- fields

  df <- data.frame(values[columns], check.names = FALSE)
  if( length(index_col) && index_col %in% names(values) ){
    rownames(df) <- make.unique(as.character(values[[index_col]]))
  }

  df
}