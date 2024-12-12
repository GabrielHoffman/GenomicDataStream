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



#' Fit series of linear regression models
#'
#' Fit series of linear regression models
#'
#' @param y response vector
#' @param design design matrix, mat or sp_mat
#' @param data \code{matrix} or \code{GenomicDataStream} with additional features to be fit one at a time
#' @param weights sample-level weights
#' @param detail return model with specified level of detail. LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
#' @param preprojection default true. Use preproject of design matrix to accelerate calculations
#' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across features
#' @param ... other args
#'
#' @examples
#' # create response, design and weights
#' y <- seq(60)
#' X_design <- matrix(1, 60, 1)
#' w <- rep(1, 60)
#'
#' # VCF file
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' # Read data into R
#' # then run lmFitFeatures()
#' gds <- GenomicDataStream(file, "DS", initialize = TRUE)
#' dat <- getNextChunk(gds)
#' X_features <- dat$X
#'
#' res1 <- lmFitFeatures(y, X_design, X_features, w)
#'
#' # Data stays at C++ level
#' # then run lmFitFeatures()
#' gds <- GenomicDataStream(file, "DS")
#'
#' res2 <- lmFitFeatures(y, X_design, gds, w)
#
#' @importMethodsFrom fastLinReg lmFitFeatures
#' @export
#' @rdname lmFitFeatures
#' @aliases lmFitFeatures,GenomicDataStream-method
setMethod(
  "lmFitFeatures", signature(data = "GenomicDataStream"),
  function(y, design, data, weights, detail = 1, preprojection = TRUE, nthreads = 1, ...) {
    if( isInitialized(data) ){
      stop("GenomicDataStream is already initialized")
    }

    if (detail > 4) stop("detail > 4 not defined")

    # if weights is not given, set to empty vector
    if (missing(weights)) {
      weights <- rep(1, 0)
    }

    # check dimensions
    stopifnot(length(y) == nrow(design))
    if (length(weights) > 0) {
      stopifnot(length(y) == length(weights))
    }

    lmFitFeatures_export(y, design, as.list(data), weights, detail, preprojection, nthreads)
  }
)











#' Fit series of linear regression models to multiple responses with shared design matrix
#'
#' Fit regression model \code{Y[j,] ~ X} for each feature j
#'
#' @param Y matrix of responses as __rows__
#' @param design design matrix
#' @param Weights matrix sample-level weights the same dimension as Y
#' @param detail level of model detail returned, with LOW = 0, MEDIUM = 1, HIGH = 2. LOW (\code{beta}, \code{se}, \code{sigSq}, \code{rdf}), MEDIUM (\code{vcov}), HIGH (\code{residuals}), MOST (\code{hatvalues})
#' @param chunkSize number of features to read per chunk
#' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
#' @param ... other args
#'
#' @details Since the weights vary for each response, each model is computed separately without recycling precomputed values
#'
#' @return List of parameter estimates with entries \code{coef},  \code{se}, \code{sigSq}, \code{rdf} and other depending on \code{detail}
#'
#' @name lmFitResponses
#' @examples
#' library(DelayedArray)
#' n <- 100
#' m <- 5
#' nc <- 2
#' set.seed(1)
#' Y <- matrix(rnorm(n * m), m, n)
#' Y <- DelayedArray(Y)
#' X <- matrix(rnorm(n * nc), n, nc)
#' rownames(Y) <- seq(m)
#' W <- matrix(runif(n * m), m, n)
#'
#' # fit regressions with model j using Y[j,] as a response
#' fit <- lmFitResponses(Y, X, W)
#' #
#' # examine results
#' lapply(fit, head, 2)
#
#' @importMethodsFrom fastLinReg lmFitResponses
#' @importFrom beachmat initializeCpp
#' @export
#' @rdname lmFitResponses
#' @aliases lmFitResponses,ANY-method
setMethod(
  "lmFitResponses", signature(Y = "ANY"),
  function(Y, design, Weights, detail = 1, chunkSize = 1000, nthreads = 1, ...) {

    if (detail > 4) stop("detail > 4 not defined")

    # if weights is not given, set to empty vector
    if (missing(Weights)) {
      Weights <- matrix(1, nrow(Y), ncol(Y))
    }

    # check dimensions
    stopifnot(ncol(Y) == nrow(design))
    stopifnot(all(dim(Y) == dim(Weights)))

    ids <- rownames(Y)
    if (is.null(ids)) {
      ids <- paste0("resp_", seq(nrow(Y)))
    }

    ptr <- initializeCpp(Y)

    lmFitResponses_export(ptr, design, ids, Weights, chunkSize, detail, nthreads, verbose = TRUE)
  }
)
