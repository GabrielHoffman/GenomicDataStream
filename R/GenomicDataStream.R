#' GenomicDataStream
#'
#' Read genomic data files (VCF, BCF, BGEN, h5ad) into R/Rcpp in chunks for analysis with Armadillo or Eigen libraries
#'
#' @name GenomicDataStream
#' @useDynLib GenomicDataStream
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
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
  MAF = "numeric", 
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
#' @param field field of VCF/BCF to read.  Ignored for other file types
#' @param region target in the format \code{chr2:1-12345}. Multiple regions can be separated by one of \code{",\n\t"}, for example \code{"chr2:1-12345, chr3:1000-8000"}. Setting region to \code{""} includes all variants
#' @param samples string of comma separated sample IDs to extract: \code{"ID1,ID2,ID3"}.  \code{"-"} indicates all samples
#' @param MAF minor allele frequency filter applied to variants with max value <= 2. Use \code{NaN} to retain all variants
#' @param minVariance minimum variance filter applied to variants with max value > 2
#' @param chunkSize	number of variants to return per chunk
#' @param missingToMean	if true, set missing values to the mean dosage value. if false, set to \code{NaN}
#' @param initialize default \code{FALSE}.  If \code{TRUE}, file info is read from path, otherwise store path until \code{GenomicDataStream} is initialized later
#'
#' @details Variants are filtered using \code{MAF} if the max value is <=2, or \code{minVariance} otherwise
#
# @details Consider minor allele frequency (MAF) \eqn{f} and  Hardy-Weinberg equilibrium, the allelic states have probability \eqn{f^2, 2f(1-f), (1-f)^2}. If the variant has mean \eqn{\mu} and variance \eqn{\sigma^2}, MAF can be estimated from the mean as \eqn{min(\mu/2, 1 - \mu/2)} or from the variance as \eqn{p = 1+sqrt(1-2\sigma^2))/2} and MAF = \eqn{min(p, 1-p)}.  In addition the sample variance of the variant is \eqn{2(1-f)f}.  Therefore, setting a MAF cutoff corresponds to a variance cutoff that can also apply to multi-allelic variants.
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
#' 
#' @importFrom methods new is
#' @export
GenomicDataStream <- function(file, field = "", region = "", samples = "-", MAF = 0, minVariance = 0, chunkSize = 10000, missingToMean = TRUE, initialize = FALSE){

  chunkSize <- as.integer(chunkSize)
  samples <- paste(samples, collapse = ",")

  file = path.expand(file)

  # check that file exists
  if( ! file.exists(file) ){
    stop("File does not exist")
  }

  if( !is.nan(MAF) & MAF < 0 ){
    stop("MAF must be >= 0")
  }
  if( !is.nan(minVariance) & minVariance < 0 ){
    stop("minVariance must be >= 0")
  }

  if (initialize) {
    # Create GenomicDataStream and return external pointer
    ptr <- create_xptr(file, field, region, samples, MAF, minVariance, chunkSize, missingToMean)

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
      MAF = MAF,
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
      MAF = MAF,
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
#' 
#' @export
isInitialized <- function(x) {
  stopifnot(is(x, "GenomicDataStream"))
  x@initialized
}

#' Initialize GenomicDataStream
#'
#' Read file info from path to initialise stream. If already initialized, return to the beginning of the stream
#'
#' @param x \code{GenomicDataStream}
#' @param region target in the format \code{chr2:1-12345}. Multiple regions can be separated by one of \code{",\n\t"}, for example \code{"chr2:1-12345, chr3:1000-8000"}. Setting region to \code{""} includes all variants
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
#' 
#' @export
initializeStream <- function(x, region = NULL) {

  if(is.null(region)) region <- x@region

  if (isInitialized(x)) {
    x <- setRegion(x, region)
    return(x)
  }

  # Create initialized GenomicDataStream
  GenomicDataStream(
    file = x@file,
    field = x@field,
    region = x@region,
    samples = x@samples,
    MAF = x@MAF,
    minVariance = x@minVariance,
    chunkSize = x@chunkSize,
    missingToMean = x@missingToMean,
    initialize = TRUE
  )
}


#' Set Chunk Size
#'
#' Set chunk size for existing \code{GenomicDataStream}
#'
#' @param x \code{GenomicDataStream}
#' @param chunkSize positive integer
#'
#' @return none
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize=TRUE)
#'
#' obj <- setChunkSize(obj, 200)
#' 
#' @export
setChunkSize <- function (x, chunkSize) {

  chunkSize <- as.integer(chunkSize)
  
  if ( isInitialized(x) ) {
    ptr <- setChunkSize_rcpp(x@ptr, chunkSize)

    # get additional information about data
    info <- getInfo(ptr)

    obj <- new("GenomicDataStream",
        initialized = TRUE,
        ptr = ptr,
        file = x@file,
        field = x@field,
        region = x@region,
        samples = x@samples,
        MAF = x@MAF,
        minVariance = x@minVariance,
        chunkSize = chunkSize,
        missingToMean = x@missingToMean,
        streamType = info$streamType,
        nsamples = info$nsamples)
  }else{
    obj <- new("GenomicDataStream",
      initialized = FALSE,
      file = x@file,
      field = x@field,
      region = x@region,
      samples = x@samples,
      MAF = x@MAF,
      minVariance = x@minVariance,
      chunkSize = chunkSize,
      missingToMean = x@missingToMean)
  }

  obj 
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
        MAF = x@MAF,
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
      MAF = x@MAF,
      minVariance = x@minVariance,
      chunkSize = x@chunkSize,
      missingToMean = x@missingToMean)
  }

  obj 
}



#' Get sample names 
#'
#' Get sample names in order that the genotypes are extracted
#'
#' @param x \code{GenomicDataStream}
#'
#' @return array of string names
#'
#' @description BGEN uses sample order from the query, but VCF/BCF/PGEN uses order in file
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' obj <- GenomicDataStream(file, "DS", init=TRUE)
#'
#' getSampleNames(obj)
#
#' @export
getSampleNames <- function(x) {
  
  if ( ! isInitialized(x) ) {
    stop("Must be initialized first")
  }

  return( getSampleNames_rcpp( x@ptr ) )
}

#' @export
#' @docType methods
#' @keywords internal
#' @aliases rownames
setGeneric("rownames", function(x, do.NULL = TRUE, prefix = "row") standardGeneric("rownames"))

#' Get rownames
#'
#' Get rownames (i.e. sample names) in order that the samples are extracted
#'
#' @param x \code{GenomicDataStream}
#' @param do.NULL not used
#' @param prefix not used
#'
#' @return array of string names
#'
#' @description BGEN uses sample order from the query, but VCF/BCF/PGEN uses order in file
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' obj <- GenomicDataStream(file, "DS", init=TRUE)
#'
#' rownames(obj)
#
#' @rawNamespace export(rownames)
#' @rdname rownames
#' @aliases rownames,GenomicDataStream-method
setMethod(
  "rownames", signature(x = "GenomicDataStream"),
  function(x, do.NULL = TRUE, prefix = "row"){
  getSampleNames(x)
})




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
featuresRead <- function(x) {
  stopifnot(is(x, "GenomicDataStream"))
  stopifnot(isInitialized(x))

  featuresRead_rcpp(x@ptr)
}

#' Get ranges for each chromosome
#' 
#' Get max and max postion for each chromosome
#' 
#' @param x \code{GenomicDataStream}
#' 
#' @examples
#' file <- system.file("extdata", "test.bed", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, chunkSize = 5, initialize = TRUE)
#' getChromRanges( obj )
#
#' @return \code{data.frame} of chrom, start, end
#' @export
getChromRanges <- function( x ){

  stopifnot(is(x, "GenomicDataStream"))
  stopifnot(isInitialized(x))

  getChromRanges_rcpp(x@ptr)
}

#' Get genomic intervals
#' 
#' Get genomic intervals chopping up regions into smaller chunks
#' 
#' @param x \code{data.frame} with columns chrom, start, end
#' @param nchunks number of chunks to create
#' 
#' @examples
#' file <- system.file("extdata", "test.bed", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, chunkSize = 5, initialize = TRUE)
#' df <- getChromRanges(obj)
#' chopChroms( df )
#
#' @export
chopChroms <- function( x, nchunks = 10){

  stopifnot(nchunks > 2)

  df <- x

  if( nrow(df) == 0 ){
    return("")
  }

  k <- (nchunks+1) / seq(nrow(df))

  regions <- sapply(seq(nrow(df)), function(i){

    pos <- seq(df$start[i], df$end[i], length.out=k)
    pos <- floor(pos)

    regions <- c()
    for(j in 2:length(pos)){
      offset <- (j!=2)
      reg <- paste0(df$chrom[i], ":", pos[j-1] + offset, "-", pos[j])
      regions <- c(regions, reg)
    }
    regions
    })
  c(regions)
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
      if( object@streamType %in% c("vcf", "vcf.gz", "bcf")){
        cat("  field:        ", object@field, "\n")
      }
      cat("  region:       ", object@region, "\n")
      cat("  samples:      ", object@nsamples, "\n")
      cat("  MAF:          ", object@MAF, "\n")
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


#' Get location of each feature
#'
#' Get location of each feature in GenomicDataStream
#'
#' @param x \code{GenomicDataStream}
#'
#' @return get data chunk as \code{list} with entries \code{X} storing a matrix with features as columns, and \code{info} storing information about each feature as rows
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)
#'
#' getVariantLocations(obj)
#' 
#' @export
getVariantLocations = function(x){

  # pass R CMD check
  CHROM <- POS <- NULL

  x <- initializeStream(x)

  lst <- list()
  i <- 1

  # loop until break
  while (1) {
    # get data chunk
    # data$X matrix with features as columns
    # data$info information about each feature as rows
    dat <- getNextChunk(x)

    if (atEndOfStream(x)) break

    lst[[i]] <- dat$info
    i <- i + 1
  }
  lst <- do.call(rbind, lst)

  # locations as BED coordinations
  with(lst, paste0(CHROM, ":", POS, "-", POS))
}


#' Get chunk sizes
#' 
#' Get array of chunk sizes
#' 
#' @param x \code{GenomicDataStream}
#' 
#' @return array of chunk sizes
#' 
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, "DS", chunkSize = 5, initialize = TRUE)
#'
#' countChunks(obj)
#' 
#' @export
countChunks = function(x){
  counts <- c()

  # count number of variants
  x <- initializeStream(x)
  while(1){
    dat <- getNextChunk(x)
    counts <- append(counts, ncol(dat$X))
    if (atEndOfStream(x)) break
  }
  counts
}


#' Collapse array of regions
#'
#' Collapse array of regions into one interval per chromosome
#' 
#' @param regions array of interval strings
#' 
#' @return \code{data.frame} one interval per chromosome spanning all given regions
#' 
#' @examples
#' regions = c( "chr2:3-5", "chr1:4-5", "chr1:1-3")
#' collapseRegions( regions )
#
#' @importFrom stringr str_split
#' @importFrom dplyr arrange tibble group_by summarize
#' @importFrom magrittr `%>%`
#' @export
#' @keywords internal
collapseRegions = function(regions){

  start <- end <- NULL

  res1 <- str_split(regions, ":")

  chrom <- sapply(res1, function(x) x[1])
  pos <- sapply(res1, function(x) x[2])
  res2 <- str_split(pos, "-")

  tibble(chrom = chrom, 
        start = as.numeric(sapply(res2, function(x) x[1])),
        end = as.numeric(sapply(res2, function(x) x[2]))) %>%
        group_by(chrom) %>%
        summarize(data.frame(chrom = chrom[1], 
                          start = min(start), 
                          end = max(end))) %>%
        as.data.frame 
}










