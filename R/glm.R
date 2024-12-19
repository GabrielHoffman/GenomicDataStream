

#' Fit series of generalized linear regression models
#'
#' Fit series of generalized linear regression models
#'
#' @param y response vector
#' @param design design matrix shared across all models
#' @param data \code{matrix} or \code{GenomicDataStream} with additional features to be fit one at a time
#' @param family a description of the error distribution and link function to be used in the modelm just like for \code{glm()}.  Also supports negative binomial as string \code{"nb:theta"}, see details below
#' @param weights vector of sample-level weights
#' @param offset vector of sample-level offset values
#' @param detail level of model detail returned, with LEAST = 0, LOW = 1, MEDIUM = 2, HIGH = 3, MOST = 4. LEAST (beta), LOW (beta, se, dispersion, rdf), MEDIUM (vcov), HIGH (pearson residuals), MOST (hatvalues)
#' @param doCoxReid use Cox-Reid adjustment when estimating overdispersion for negative binomial models.  Default TRUE for less than 100 samples
#' @param shareTheta estimate theta from design matrix, and share across all features instead of re-estimating for each feature
#' @param fastApprox default false.  if true, use pre-projection on the working response from an initial regression fit on only the design.  Under the null for data, this is a very good approximation and _much_ faster
#' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across features
#' @param epsilon tolerance for GLM IRLS 
#' @param maxit max iterations for GLM IRLS
#' @param epsilon_nb tolerance for negative binomial
#' @param maxit_nb max iterations for negative binomial
#' @param ... other args
#'
#' @details 
#' Generalized linear models can be fit with \code{family} like in \code{glm()} using \code{gaussian()}, \code{poisson()}, \code{binomial()}, \code{binomial("probit")}, \code{quasibinomial()}, \code{quasipoisson()}, \code{negative.binomial(theta)}, \code{"nb"}, \code{"nb:theta"}. Or array of entries of form \code{"nb:theta"}, where \code{theta} is the parameter for the negative binomial distribution
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
#' @importFrom fastLinReg getFamilyString
#' @export
#' @rdname glmFitFeatures
#' @aliases glmFitFeatures,ANY,ANY,GenomicDataStream-method
setMethod(
  "glmFitFeatures", signature(data = "GenomicDataStream"),
  function(y, design, data, family, weights, offset, detail = 1, doCoxReid = length(y) < 100, shareTheta = FALSE, fastApprox = FALSE, nthreads = 1, epsilon = 1e-8, maxit = 25, epsilon_nb = 1e-4, maxit_nb = 5, ...) {
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

    family = getFamilyString(family)

    binom = c("binomial/logit", "binomial/probit")

    if( fastApprox && family %in% binom && ! all(y %in% c(0,1)) ){
      stop("fastApprox performs poorly with binomial family for fractional response")
    }

    glmFitFeatures_export(y, design, as.list(data), family, weights, offset, detail, doCoxReid, shareTheta,fastApprox, nthreads, epsilon, maxit, epsilon_nb, maxit_nb,...)
  }
)











#' Fit series of generalized linear models to multiple responses with shared design matrix
#'
#' Fit regression model \code{Y[j,] ~ X} for each feature j
#'
#' @param Y matrix of responses as __rows__
#' @param design design matrix
#' @param family a description of the error distribution and link function to be used in the modelm just like for \code{glm()}.  Also supports negative binomial as string \code{"nb:theta"}, see details below 
#' @param weights vector of sample-level weights
#' @param offset vector of sample-level offset values
#' @param detail level of model detail returned, with LEAST = 0, LOW = 1, MEDIUM = 2, HIGH = 3, MOST = 4. LEAST (beta), LOW (beta, se, dispersion, rdf), MEDIUM (vcov), HIGH (pearson residuals), MOST (hatvalues)
#' @param doCoxReid use Cox-Reid adjustment when estimating overdispersion for negative binomial models.  Default TRUE for less than 100 samples
#' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
#' @param epsilon tolerance for GLM IRLS 
#' @param maxit max iterations for GLM IRLS
#' @param epsilon_nb tolerance for negative binomial
#' @param maxit_nb max iterations for negative binomial
#' @param chunkSize number of features to read per chunk
#' @param ... other args
#'
#' @details 
#' Generalized linear models can be fit with \code{family} like in \code{glm()} using \code{gaussian()}, \code{poisson()}, \code{binomial()}, \code{binomial("probit")}, \code{quasibinomial()}, \code{quasipoisson()}, \code{negative.binomial(theta)}, \code{"nb"}, \code{"nb:theta"}. Or array of entries of form \code{"nb:theta"}, where \code{theta} is the parameter for the negative binomial distribution
#'
#' @return List of parameter estimates with entries \code{coef},  \code{se}, \code{dispersion}, \code{rdf} and other depending on \code{detail}
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
#' @importFrom fastLinReg getFamilyString
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

    family = getFamilyString( family )

    # if family is a string,
    # repeat for each response
    if( length(family) == 1){
      family = rep(family, nrow(Y))
    }else{
      # entries in family are only allowed to vary
      # if nb:numeric is used
      res = grepl("^nb:(\\d+)", family)

      if( length(unique(family)) != 1 && any(!res) ){
        stop("family can only vary if nb:xx is used")
      } 
    }

    if( length(family) != nrow(Y)){
      stop("family must be either a string or an array with one entry per response")
    }

    ptr <- initializeCpp(Y)

    glmFitResponses_export(ptr, design, ids, family, weights, offset, chunkSize, detail, doCoxReid, nthreads, epsilon, maxit, epsilon_nb, maxit_nb)
  }
)
