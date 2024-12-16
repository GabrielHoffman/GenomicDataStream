

#' Fit series of generalized linear regression models
#'
#' Fit series of generalized linear regression models
#'
#' @param y response vector
#' @param design design matrix shared across all models
#' @param data \code{matrix} or \code{GenomicDataStream} with additional features to be fit one at a time
#' @param family family type of GLM: "gaussian", "logit", "probit", "poisson", "nb:x" where x is a numeric value of theta 
#' @param weights vector of sample-level weights
#' @param offset vector of sample-level offset values
#' @param detail level of model detail returned, with LEAST = 0, LOW = 1, MEDIUM = 2, HIGH = 3, MOST = 4. LEAST (beta), LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
#' @param doCoxReid use Cox-Reid adjustment when estimating overdispersion for negative binomial models.  Default TRUE for less than 100 samples
#' @param shareTheta estimate theta from design matrix, and share across all features instead of re-estimating for each feature
#' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across features
#' @param epsilon tolerance for GLM IRLS 
#' @param maxit max iterations for GLM IRLS
#' @param epsilon_nb tolerance for negative binomial
#' @param maxit_nb max iterations for negative binomial
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
#' # then run glmFitFeatures()
#' gds <- GenomicDataStream(file, "DS", initialize = TRUE)
#' dat <- getNextChunk(gds)
#' X_features <- dat$X
#'
#' res1 <- glmFitFeatures(y, X_design, X_features, family="gaussian", w)
#'
#' # Data stays at C++ level
#' # then run lmFitFeatures()
#' gds <- GenomicDataStream(file, "DS")
#'
#' res2 <- glmFitFeatures(y, X_design, gds, family="gaussian", w)
#
#' @importMethodsFrom fastLinReg glmFitFeatures
#' @export
#' @rdname glmFitFeatures
#' @aliases glmFitFeatures,ANY,ANY,GenomicDataStream-method
setMethod(
  "glmFitFeatures", signature(data = "GenomicDataStream"),
  function(y, design, data, family, weights, offset, detail = 1, doCoxReid = length(y) < 100, shareTheta = FALSE, nthreads = 1, epsilon = 1e-8, maxit = 25, epsilon_nb = 1e-4, maxit_nb = 5, ...) {
    if( isInitialized(data) ){
      stop("GenomicDataStream is already initialized")
    }

    if (detail > 4) stop("detail > 4 not defined")

   # if weights is not given, set to empty vector
    if (missing(weights)) {
      weights <- rep(1, 0)
    }
    if (missing(offset)) {
      offset <- rep(0, length(y))
    }

    # check dimensions
    stopifnot(length(y) == nrow(design))
    if (length(weights) > 0) {
      stopifnot(length(y) == length(weights))
    }

    glmFitFeatures_export(y, design, as.list(data), family, weights, offset, detail, doCoxReid, shareTheta, nthreads, epsilon, maxit, epsilon_nb, maxit_nb,...)
  }
)











#' Fit series of generalized linear models to multiple responses with shared design matrix
#'
#' Fit regression model \code{Y[j,] ~ X} for each feature j
#'
#' @param Y matrix of responses as __rows__
#' @param design design matrix
#' @param family type of GLM: "gaussian", "logit", "probit", "poisson", "nb:x" where x is a numeric value of theta  
#' @param weights vector of sample-level weights
#' @param offset vector of sample-level offset values
#' @param detail level of model detail returned, with LEAST = 0, LOW = 1, MEDIUM = 2, HIGH = 3, MOST = 4. LEAST (beta), LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
#' @param doCoxReid use Cox-Reid adjustment when estimating overdispersion for negative binomial models.  Default TRUE for less than 100 samples
#' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
#' @param epsilon tolerance for GLM IRLS 
#' @param maxit max iterations for GLM IRLS
#' @param epsilon_nb tolerance for negative binomial
#' @param maxit_nb max iterations for negative binomial
#' @param chunkSize number of features to read per chunk
#' @param ... other args
#'
#' @return List of parameter estimates with entries \code{coef},  \code{se}, \code{sigSq}, \code{rdf} and other depending on \code{detail}
#'
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
#'
#' # fit regressions with model j using Y[j,] as a response
#' fit <- glmFitResponses(Y, X, "gaussian")
#' 
#' # examine results
#' lapply(fit, head, 2)
#
#' @importMethodsFrom fastLinReg glmFitResponses
#' @importFrom beachmat initializeCpp
#' @export
#' @rdname glmFitResponses
#' @aliases glmFitResponses,ANY-method
setMethod(
  "glmFitResponses", signature(Y = "ANY"),
  function(Y, design, family, weights, offset, detail = 1, doCoxReid = nrow(design) < 100, nthreads = 1, epsilon = 1e-8, maxit = 25, epsilon_nb = 1e-4, maxit_nb = 5, chunkSize = 1000, ...) {
   
    if (detail > 4) stop("detail > 4 not defined")

    # if weights is not given, set to empty vector
    if (missing(weights)) {
      weights <- rep(1, nrow(design))
    }
    if (missing(offset)) {
      offset <- rep(0, nrow(design))
    }

    # check dimensions
    stopifnot(ncol(Y) == nrow(design))

    ids <- rownames(Y)
    if (is.null(ids)) {
      ids <- paste0("resp_", seq(nrow(Y)))
    }

    ptr <- initializeCpp(Y)

    glmFitResponses_export(ptr, design, ids, family, weights, offset, chunkSize, detail, doCoxReid, nthreads, epsilon, maxit, epsilon_nb, maxit_nb)
  }
)
