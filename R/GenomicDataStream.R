
#' GenomicDataStream
#'
#' Read genomic data files (VCF, BCF, BGEN, h5ad) into R/Rcpp in chunks for analysis with Armadillo or Eigen libraries
#'
#' @name GenomicDataStream
#' @useDynLib GenomicDataStream 
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
NULL


setClass("GenomicDataStream", slots=list(ptr="externalptr", file="character", field = "character", region="character", samples="character", chunkSize="integer", missingToMean="logical", featuresRead="integer", streamType="character", nsamples="integer"))

#' @export
GenomicDataStream = function( file, field = "", region = "", samples="-", chunkSize = 1000, missingToMean = FALSE){

	chunkSize = as.integer(chunkSize)
	samples = paste(samples, collapse=',')

	# Create GenomicDataStream and return external pointer
	ptr = create_xptr( file, field, region, samples, chunkSize, missingToMean)

	# get additional information about data
	info = getInfo( ptr )

	# return object
	new( "GenomicDataStream", ptr = ptr, 
								file = file,
								field = field,
								region = region,
								samples = samples,
								chunkSize = chunkSize,
								missingToMean = missingToMean,
								featuresRead = 0L,
								streamType = info$streamType,
								nsamples = info$nsamples)
}


#' @export
hasReachedEnd = function( x ){

	stopifnot( is(x, "GenomicDataStream") )

	hasReachedEnd_rcpp(x@ptr)
}

#' @export
featuresRead = function( x ){

	stopifnot( is(x, "GenomicDataStream") )

	featuresRead_rcpp(x@ptr)
}

#' @export
getNextChunk = function( x ){

	stopifnot( is(x, "GenomicDataStream") )

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
#' @importFrom S4Vectors coolcat
#' @aliases show,GenomicDataStream,GenomicDataStream-method
#' @export
setMethod(
  "show", "GenomicDataStream",
  function(object) {
    cat("\t\t", class(object), "\n\n")
    cat("  file:", basename(object@file), "\n")
    cat("  stream type:", object@streamType, "\n")
    cat("  field:", object@field, "\n")
    cat("  region:", object@region, "\n")

    # if( object@samples == "-"){
    # 	cat("samples: -\n")
    # }else{
    # 	ids = strsplit(object@samples, ',')[[1]]
	#     coolcat("samples(%d): %s\n", ids)
	# }
	cat("  samples:", obj@nsamples, "\n") 

    cat("  missingToMean:", object@missingToMean, "\n") 
    cat("  chunkSize:", object@chunkSize, "\n")
    cat("  features read:", object@featuresRead, "\n") 
    cat("  reached end:", hasReachedEnd(object), "\n")  
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
setMethod("print", "GenomicDataStream", function(x,...)
	show(x)
	)