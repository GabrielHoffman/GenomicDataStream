#' GenomicDataStream
#'
#' Read genomic data files (VCF, BCF, BGEN, h5ad) into R/Rcpp in chunks for analysis with Armadillo or Eigen libraries
#'
#' @name GenomicDataStream
#' @useDynLib GenomicDataStream
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
NULL

#' Interface to genomic data files
#'
#' Read genomic data files (VCF, BCF, BGEN, h5ad) into R/Rcpp in chunks for analysis with Armadillo or Eigen libraries
#'
#' @export
setClass("GenomicDataStream", slots = list(
  initialized = "logical", 
  ptr = "externalptr", 
  file = "character", 
  field = "character", 
  region = "character", 
  samples = "character", 
  chunkSize = "integer", 
  missingToMean = "logical", 
  featuresRead = "integer", 
  streamType = "character", 
  nsamples = "integer", 
  minVariance = "numeric"))


#' @export
as.list.GenomicDataStream <- function(x, ...) {
  list(
    file = x@file,
    field = x@field,
    region = x@region,
    samples = x@samples,
    chunkSize = x@chunkSize,
    missingToMean = x@missingToMean
  )
}

#' GenomicDataStream to read genotypes / dosages
#'
#' Interface to GenomicDataStream C++ code
#'
#' @param file file in VCF/BCF/BGEN/PGEN format with index
#' @param field field of VCF/BCF to read
#' @param region target in the format \code{chr2:1-12345}. Multiple regions can be separated by one of \code{",\n\t"}, for example \code{"chr2:1-12345, chr3:1000-8000"}. Setting region to \code{""} includes all variants
#' @param samples string of comma separated sample IDs to extract: \code{"ID1,ID2,ID3"}.  \code{"-"} indicates all samples
#' @param MAF generalized minor allele frequency cutoff keeps variants with variance > \eqn{2(1-f)f}.
#' @param minVariance features with variance \code{>= minVariance} are retained.  Defaults to \eqn{2(1-f)f}. If \code{NaN}, no filtering is applied.
#' @param chunkSize	number of variants to return per chunk
#' @param missingToMean	if true, set missing values to the mean dosage value. if false, set to \code{NaN}
#' @param initialize default \code{FALSE}.  If \code{TRUE}, file info is read from path, otherwise store path until \code{GenomicDataStream} is initialized later
#'
#' @details Consider minor allele frequency (MAF) \eqn{f} and  Hardy-Weinberg equilibrium, the allelic states have probability \eqn{f^2, 2f(1-f), (1-f)^2}. If the variant has mean \eqn{\mu} and variance \eqn{\sigma^2}, MAF can be estimated from the mean as \eqn{min(\mu/2, 1 - \mu/2)} or from the variance as \eqn{p = 1+sqrt(1-2\sigma^2))/2} and MAF = \eqn{min(p, 1-p)}.  In addition the sample variance of the variant is \eqn{2(1-f)f}.  Therefore, setting a MAF cutoff corresponds to a variance cutoff that can also apply to multi-allelic variants.
#'
#' @return object of class \code{GenomicDataStream}
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)
#'
#' obj
#'
#' # loop until break
#' while (1) {
#'   # get data chunk
#'   # data$X matrix with features as columns
#'   # data$info information about each feature as rows
#'   dat <- getNextChunk(obj)
#'
#'   if (atEndOfStream(obj)) break
#'
#'   print(dat$info)
#' }
#' #
#' @importFrom methods new is
#' @export
GenomicDataStream <- function(file, field = "", region = "", samples = "-", MAF = 0, minVariance = 2*(1-MAF)*MAF, chunkSize = 10000, missingToMean = TRUE, initialize = FALSE){

  chunkSize <- as.integer(chunkSize)
  samples <- paste(samples, collapse = ",")

  file = path.expand(file)

  # check that file exists
  if( ! file.exists(file) ){
    stop("File does not exist")
  }

  if( MAF < 0 ){
    stop("MAF must be >= 0")
  }
  if( !is.nan(minVariance) & minVariance < 0 ){
    stop("minVariance must be >= 0")
  }

  if (initialize) {
    # Create GenomicDataStream and return external pointer
    ptr <- create_xptr(file, field, region, samples, minVariance, chunkSize, missingToMean)

    # get additional information about data
    info <- getInfo(ptr)

    # return object
    obj <- new("GenomicDataStream",
      initialized = TRUE,
      ptr = ptr,
      file = file,
      field = field,
      region = region,
      samples = samples,
      minVariance = minVariance,
      chunkSize = chunkSize,
      missingToMean = missingToMean,
      streamType = info$streamType,
      nsamples = info$nsamples
    )
  } else {
    # return object
    obj <- new("GenomicDataStream",
      initialized = FALSE,
      file = file,
      field = field,
      region = region,
      samples = samples,
      minVariance = minVariance,
      chunkSize = chunkSize,
      missingToMean = missingToMean
    )
  }

  obj
}

#' Get status of GenomicDataStream
#'
#' If \code{initialized}, return \code{TRUE}, else \code{FALSE}
#'
#' @param x \code{GenomicDataStream}
#'
#' @return initialization status
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5)
#'
#' # by default, GenomicDataStream is not initialized
#' isInitialized(obj)
#'
#' # initialize
#' obj <- initializeStream(obj)
#'
#' isInitialized(obj)
#' #
#' @export
isInitialized <- function(x) {
  stopifnot(is(x, "GenomicDataStream"))
  x@initialized
}

#' Initialize GenomicDataStream
#'
#' Read file info from path to initialise stream
#'
#' @param x \code{GenomicDataStream}
#'
#' @return initialized \code{GenomicDataStream}
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5)
#'
#' # by default, GenomicDataStream is not initialized
#' isInitialized(obj)
#'
#' # initialize
#' obj <- initializeStream(obj)
#'
#' isInitialized(obj)
#' #
#' @export
initializeStream <- function(x) {
  if (isInitialized(x)) {
    return(x)
  }

  # Create initialized GenomicDataStream
  GenomicDataStream(
    file = x@file,
    field = x@field,
    region = x@region,
    samples = x@samples,
    minVariance = x@minVariance,
    chunkSize = x@chunkSize,
    missingToMean = x@missingToMean,
    initialize = TRUE
  )
}


#' Set regions of GenomicDataStream
#'
#' Set regions of GenomicDataStream
#'
#' @param x \code{GenomicDataStream}
#' @param region target in the format \code{chr2:1-12345}. Multiple regions can be separated by one of \code{",\n\t"}, for example \code{"chr2:1-12345, chr3:1000-8000"}. Setting region to \code{""} includes all variants
#'
#' @return \code{GenomicDataStream} with region set
#'
#' @description If \code{GenomicDataStream} is already initialized, set new query in C++ backend. Otherwise substitute \code{region} values
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5)
#'
#' # by default, GenomicDataStream is not initialized
#' setRegion(obj, "1:10000-12000")
#
#' @export
setRegion <- function(x, region) {
  
  if ( isInitialized(x) ) {
    ptr <- setRegions_rcpp(x@ptr, region)

    # get additional information about data
    info <- getInfo(ptr)

    obj <- new("GenomicDataStream",
        initialized = TRUE,
        ptr = ptr,
        file = x@file,
        field = x@field,
        region = region,
        samples = x@samples,
        minVariance = x@minVariance,
        chunkSize = x@chunkSize,
        missingToMean = x@missingToMean,
        streamType = info$streamType,
        nsamples = info$nsamples)
  }else{
    obj <- new("GenomicDataStream",
      initialized = FALSE,
      file = x@file,
      field = x@field,
      region = region,
      samples = x@samples,
      minVariance = x@minVariance,
      chunkSize = x@chunkSize,
      missingToMean = x@missingToMean)
  }

  obj 
}

#' Detected if end of stream is reached
#'
#' Detected if end of stream is reached
#'
#' @param x \code{GenomicDataStream}
#'
#' @return if end of stream has been reached, return \code{TRUE}.  Else \code{FALSE}
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)
#'
#' # loop until break
#' while (1) {
#'   # get data chunk
#'   # data$X matrix with features as columns
#'   # data$info information about each feature as rows
#'   dat <- getNextChunk(obj)
#'
#'   if (atEndOfStream(obj)) break
#'
#'   print(dat$info)
#' }
#' @export
atEndOfStream <- function(x) {
  stopifnot(is(x, "GenomicDataStream"))
  stopifnot(isInitialized(x))

  atEndOfStream_rcpp(x@ptr)
}

#' Get number of features read from GenomicDataStream
#'
#' Get number of total features read from GenomicDataStream
#'
#' @param x \code{GenomicDataStream}
#'
#' @return  total number of features read from the stream
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)
#'
#' # loop until break
#' while (1) {
#'   # get data chunk
#'   # data$X matrix with features as columns
#'   # data$info information about each feature as rows
#'   dat <- getNextChunk(obj)
#'
#'   if (atEndOfStream(obj)) break
#'
#'   print(dat$info)
#' }
#'
#' featuresRead(obj)
#' @export
#' @export
featuresRead <- function(x) {
  stopifnot(is(x, "GenomicDataStream"))
  stopifnot(isInitialized(x))

  featuresRead_rcpp(x@ptr)
}

#' Get data chunk from GenomicDataStream
#'
#' Get data chunk from GenomicDataStream
#'
#' @param x \code{GenomicDataStream}
#'
#' @return  get data chunk as \code{list} with entries \code{X} storing a matrix with features as columns, and \code{info} storing information about each feature as rows
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)
#'
#' # loop until break
#' while (1) {
#'   # get data chunk
#'   # data$X matrix with features as columns
#'   # data$info information about each feature as rows
#'   dat <- getNextChunk(obj)
#'
#'   if (atEndOfStream(obj)) break
#'
#'   print(dat$info)
#' }
#' @export
getNextChunk <- function(x) {
  stopifnot(is(x, "GenomicDataStream"))
  stopifnot(isInitialized(x))

  getNextChunk_rcpp(x@ptr)
}



#' Show object
#'
#' Show object
#'
#' @param object \code{GenomicDataStream} object
#'
#' @rdname show-methods
#' @importFrom utils head tail
#' @aliases show,GenomicDataStream,GenomicDataStream-method
#' @importFrom methods show
#' @export
setMethod(
  "show", "GenomicDataStream",
  function(object) {
    cat("\t\t", class(object), "\n\n")
    cat("  file:         ", basename(object@file), "\n")
    cat("  initialized:  ", isInitialized(object), "\n")
    if (isInitialized(object)) {
      cat("  stream type:  ", object@streamType, "\n")
      cat("  field:        ", object@field, "\n")
      cat("  region:       ", object@region, "\n")
      cat("  samples:      ", object@nsamples, "\n")
      cat("  minVar cutoff:", object@minVariance, "\n")
      cat("  missingToMean:", object@missingToMean, "\n")
      cat("  chunkSize:    ", object@chunkSize, "\n")
      cat("  features read:", featuresRead(object), "\n")
      cat("  end of stream:", atEndOfStream(object), "\n")
    }
  }
)


#' Print object
#'
#' Print object
#'
#' @param x \code{GenomicDataStream} object
#' @param ... other arguments
#'
#' @export
#' @rdname print-methods
#' @aliases print,GenomicDataStream,GenomicDataStream-method
setMethod("print", "GenomicDataStream", function(x, ...) {
  show(x)
})
