
#' @export
as.list.GenomicDataStream = function(x,...){

	list(file = x@file, 
			field = x@field,
		 	region = x@region,
		 	samples = x@samples, 
			chunkSize = x@chunkSize, 
			missingToMean = x@missingToMean)
}




# Restate setGeneric here from fastlmm
setGeneric(
  "lmFitFeatures",
  function(y, X_design, data, weights, detail = 0, preprojection = TRUE, nthreads = 1, ...) {
    standardGeneric("lmFitFeatures")
  }
)

setGeneric(
  "lmFitResponses",
  function(Y, X, ids, Weights, detail = 0, nthreads = 1, ...) {
    standardGeneric("lmFitResponses")
  }
)


#' @import fastlmm
#' @export
#' @rdname lmFitFeatures
#' @aliases lmFitFeatures, matrix-method
setMethod(
  "lmFitFeatures", signature(data = "GenomicDataStream"),
  function(y, X_design, data, weights, detail = 0, preprojection = TRUE, nthreads = 1, ...) {

	stopifnot( ! isInitialized(gds) )

	lmFitFeatures_export( y, X_design, as.list(data), weights, detail, preprojection, nthreads)
})




#' @import fastlmm
#' @export
#' @rdname lmFitResponses
#' @aliases lmFitResponses, matrix-method
setMethod(
  "lmFitResponses", signature(Y = "DelayedArray"),
  function(Y, X, ids, Weights, detail = 0, nthreads = 1, ...) {

	lmFitResponses_export( Y, X, ids, Weights, detail, nthreads)
})





