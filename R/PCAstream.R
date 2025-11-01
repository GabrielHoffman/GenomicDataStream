#' Window-based Randomized SVD 
#' 
#' @param x  \code{GenomicDataStream}, or matrix as \code{matrix}, \code{DelayedArray}. 
#'
#' @param k       integer; \cr
#'                specifies the target rank of the low-rank decomposition. \eqn{k} should satisfy \eqn{k << min(m,n)}.
#'
#' @param ... other argument to control streaming
#'
#' @param p       integer, optional; \cr
#'                number of additional power iterations (by default \eqn{p=7}).
#'
#' @param s       integer, optional; \cr
#'                oversampling parameter (by default \eqn{s=20}).
#' 
#' @param B       integer, optional; \cr
#'                number of windows (by default \eqn{B=64}).
#' 
#' @param threads integer, optional; \cr
#'                number of threads (by default \eqn{threads=4}) to read data
#'
#' @param threads2 integer, optional; \cr
#'                number of threads (by default \eqn{threads=4}), used for linear algebra opertions
#' 
#' @param scaleAndCenter bool, optional; \cr
#'                if \code{TRUE}, scale and center features
#'
#' @param shuffle  bool, optional; \cr
#'                  if \code{TRUE} (default) shuffle genomic regions, the next chunk is not in LD with the previous chunk
#' 
#' @param verbose  string, optional; \cr
#'                  if \code{TRUE} (default), print details
#' 
#' @return \code{PCAstream} returns a list containing the following three components:
#'\describe{
#'\item{d}{  array_like; \cr
#'           singular values; vector of length \eqn{(k)}.
#'}
#'
#'\item{u}{  array_like; \cr
#'           left singular vectors; \eqn{(m, k)} or \eqn{(m, nu)} dimensional array.
#'}
#'
#'\item{v}{  array_like; \cr
#'           right singular vectors; \eqn{(n, k)} or \eqn{(n, nv)} dimensional array. \cr
#'}
#'}
#'
#' @details PCAstream implements the window-based Randomized SVD proposed by Li, et al. (2023).
#' 
#' If \code{scaleAndCenter}, the data matrix is scaled so that each feature has a mean of zero and a cross-product of 1. In this case the sum of squares of all eigen values equals the number of features (i.e \code{sum(d^2) = p}).   
#'
#' Computational time is spent on two steps: 
#' 
#' 1) Reading data from disk and and processing data.  Multiple chunks can be read and processed in parallel.  This is conrolled by setting \code{threads}  
#'
#' 2) Updating PCA with current data chunk. Only one chunk can be processed at a time, but linear algebra operations can be parallelized.  This is conrolled by setting \code{threads2}  
#'
#' If \code{x} is a matrix, PCA is performed using rows as features. Note: this is the _transpose_ of what \code{svd()} does.
#'
#' @note The singular vectors are not unique and only defined up to sign. If a left singular vector has its sign changed, changing the sign of the corresponding right vector gives an equivalent decomposition.
#'
#' @references
#' \itemize{
#'  \item Li, Z., Meisner, J., & Albrechtsen, A. (2023). Fast and accurate out-of-core PCA framework for large scale biobank data. Genome Research, 33(9), 1599-1608. \doi{10.1101/gr.277525.122}.
#' }
#'
#' @author Zilong Li \email{zilong.dk@gmail.com}, Gabriel Hoffman 
#'
#' @examples
#' # Example: PCA on genotype data
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#' 
#' obj <- GenomicDataStream(file, "DS", chunkSize = 3)
#' 
#' res <- PCAstream(obj, k=5, threads=1)
#' 
#' res
#' 
#' par(pty="s")
#' plot(res)
#' 
#' @importFrom beachmat.hdf5 initializeCpp
#' @export
setGeneric(
  "PCAstream",
  function(x, k,..., p = 7, s = 20, B = 64, threads = 4, threads2 = 4, scaleAndCenter = TRUE, shuffle = TRUE, verbose = TRUE) {
    standardGeneric("PCAstream")
  }
)

#'
#' @export
#' @importFrom methods slot
#' @rdname PCAstream
#' @aliases PCAstream,GenomicDataStream-method
setMethod(
  "PCAstream", signature(x = "GenomicDataStream"),
  function(x, k,..., p = 7, s = 20, B = 64, threads = 4, threads2 = 1, scaleAndCenter = TRUE, shuffle = TRUE, verbose = TRUE) {

  stopifnot(is(x, "GenomicDataStream"))
  stopifnot(k >= 2)

  if( verbose ) message("Read through...")

  # get summary of GenomicDataStream
  sObj <- summarizeChunks(x, threads)

  N <- slot(x, "nsamples")
  M <- sum(sObj$chunks)
  s <- min(s, M-k-1)

  if( M < 2 ){
    txt <- paste("Cannot perform PCA with", M, "features")
    stop(txt)
  }

  # get query regions
  # 
  nVariants <- length(sObj$variantIDs)
  allVariantsKept <- (sObj$nVariantsBeforeFilter == nVariants)
  regions <- unique(sObj$regions)

  # if no variants were removed by filter
  if( allVariantsKept ){
    if( nVariants > 1000 ){
      # collase array of regions into smaller set of intervals
      regMerge <- collapseRegions( regions )
      regions <- chopChroms( regMerge, B)
    }
  }else{
    if( sObj$streamType %in% c("vcf.gz", "bcf") ){
      warning("Performing variant filtering and PCA on VCF/BCF files\nis substantially slower than with other formats", immediate.=TRUE)
    }
  }

  if( shuffle ){
    regions <- sample(regions, length(regions))
  }

  if( is.null(regions) ){
    stop("Chunks not read from index")
  }

  # B must be a even
  stopifnot(B %% 2 == 0)

  # set valid B
  while( ! (M > B^2) ){
    B <- B / 2
  }

  x <- setChunkSize(x, ceiling( M / B ))
  nchunks <- B

  # number of parallel chunks
  n_pll_chunks <- max(1, log2(B))
  n_pll_chunks <- min(n_pll_chunks, threads)

  # k must be < min(N,M)
  k <- min(c(k,N,M))

  if( verbose ){   
    if( length(N) > 0){ 
      message(" # samples: ", format(N, big.mark=','))
    }
    message(" # features: ", format(M, big.mark=','))
    message(" # chunks: ", format(nchunks, big.mark=','))
    message(" # PCs: ", k)
    message(" # threads: ", n_pll_chunks)
  }

  # p <- max(c(p, log2(B)+1)) 

  x <- initializeStream(x)

  # run PCA on GenomicDataStream
  res <- stream_pcaone(x@ptr, 
          region = paste(regions, collapse = "," ), 
          m = M, 
          k = k, 
          nchunks = nchunks,
          s = s, 
          p = p, 
          B = B, 
          threads = n_pll_chunks, 
          threads_eigen = threads2,
          verbose = verbose, 
          scaleAndCenter = scaleAndCenter)

  # set row and column names
  # Since variant order can be shuffled
  # res$VariantIds contains the order seen by PCA
  rownames(res$u) <- sObj$sampleIDs
  rownames(res$v) <- res$featureIds

  # if( length(res$featureIds) != length(sObj$variantIDs)){
  #   warning("Length does not match")
  # }
  res$v <- res$v[sObj$variantIDs,,drop=FALSE]
  res$featureIds <- NULL

  colnames(res$u) <- paste0("PC", seq(k))
  colnames(res$v) <- paste0("PC", seq(k))
  names(res$d) <- paste0("PC", seq(k))
  res$p <- M
  res$n <- N

  new("PCA", res)
})


#' @param chunkSize number of features to read per chunk in \code{GenomicDataStream}
#' 
#' @export
#' @importFrom methods slot
#' @importFrom beachmat initializeCpp
#' @importFrom Matrix t
#' @rdname PCAstream
#' @aliases PCAstream,ANY-method
setMethod(
  "PCAstream", signature(x = "ANY"),
  function(x, k, chunkSize = 1000, p = 7, s = 20, B = 64, threads = 4, threads2 = 4, scaleAndCenter = TRUE, shuffle = TRUE, verbose = TRUE) {

  if( scaleAndCenter ){
    if( verbose ){
      message("Scale and centering...")
    }
    # compute mean and scale of columns
    ptr <- initializeCpp(x)
    adj <- compute_center_and_scale(ptr, threads)

    # standardize columns with mean, scale and nrows
    x <- t((t(x) - adj$center) / (adj$scale*sqrt(nrow(x)-1)))
  }

  M <- nrow(x)
  N <- ncol(x)

  # k must be < min(N,M)
  k <- min(c(k,N,M))

   # B must be a even
  stopifnot(B %% 2 == 0)

  # set valid B
  while( ! (M > B^2) ){
    B <- B / 2
  }

  p <- max(c(p, log2(B)+1))   

  if( is.null(rownames(x)) ){
    rownames(x) <- paste0("f_", seq(nrow(x)))
  }  
  if( is.null(colnames(x)) ){
    colnames(x) <- paste0("s_", seq(ncol(x)))
  }

  ptr <- initializeCpp(x)

  res <- stream_pcaone_robj(
    x = ptr, 
    ids = rownames(x),
    n = N,
    chunkSize = chunkSize,
    nchunks = ceiling(M/chunkSize), 
    m = M, 
    k = k, 
    s = s, 
    p = p, 
    B = B, 
    threads = threads, 
    threads_eigen = threads2,
    verbose = verbose,
    scaleAndCenter = FALSE)

  # set row and column names
  rownames(res$u) <- colnames(x)
  rownames(res$v) <- rownames(x)

  colnames(res$u) <- paste0("PC", seq(k))
  colnames(res$v) <- paste0("PC", seq(k))
  names(res$d) <- paste0("PC", seq(k))

  res2 <- list(d = res$d, u = res$v, v = res$u)
  res2$p <- M
  res2$n <- N

  new("PCA", res2)
})


#' @param chunkSize number of features to read per chunk in \code{GenomicDataStream}
#' @param assay for \code{SummarizedExperiment} which \code{assay} to perform PCA on
#' 
#' @examples
#' # Example: PCA on gene expression data
#' library(muscat)
#' library(SingleCellExperiment)
#' library(GenomicDataStream) 
#' 
#' data(example_sce)
#' sce <- example_sce
#' 
#' # Normalize expression with log2 counts per million
#' prior.count <- 1
#' lib.size <- colSums2(counts(sce))
#' logcounts(sce) <- t(log2(t(counts(sce) + prior.count)) - log2(lib.size) + log2(1e6))
#' 
#' # PCA on SingleCellExperiment with assay logcounts
#' res <- PCAstream(sce, k=5, assay="logcounts")
#' 
#' res
#' 
#' @export
#' @importFrom methods slot
#' @importFrom beachmat initializeCpp
#' @importFrom SummarizedExperiment assay
#' @importFrom Matrix t
#' @rdname PCAstream
#' @aliases PCAstream,SummarizedExperiment-method
setMethod(
  "PCAstream", signature(x = "SummarizedExperiment"),
  function(x, k, chunkSize = 1000, assay="logcounts",..., p = 7, s = 20, B = 64, threads = 4, threads2 = 4, scaleAndCenter = TRUE, shuffle = TRUE, verbose = TRUE) {
    
  PCAstream(x = t(assay(x, assay)), 
            k = k, 
            chunkSize = chunkSize, 
            p = p, 
            s = s,
            B = B, 
            threads = threads, 
            scaleAndCenter = scaleAndCenter, 
            verbose = verbose, 
            ...)
})




#' Get information about chunks
#' 
#' Read through stream to get size of each chunk, IDs of variants and samples
#' 
#' @param x \code{GenomicDataStream}
#' @param threads number of threads
#' 
#' @return list of info about stream
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, "DS", chunkSize = 3, initialize = TRUE)
#' summarizeChunks(obj)
#' 
#' @export
summarizeChunks <- function(x, threads=4){

  x <- initializeStream(x)
  res <- summarizeChunks_rcpp(x@ptr, threads)
  res$streamType <- slot(x, "streamType")
  res
}

#' PCA result
#'
#' PCA result
#'
#' @export
setClass("PCA", contains="list")

#' Show object
#'
#' Show object
#'
#' @param object \code{PCA} object
#'
#' @return printed strings
#' @rdname show-methods
#' @importFrom utils head
#' @aliases show,PCA,PCA-method
#' @importFrom methods show
#' @export
setMethod(
  "show", "PCA",
  function(object){ 

    cat("\n       PCA: Computed first", length(object$d), "PCs\n\n")

    k <- min(3, length(object$d))

    # cat("Samples:\n")
    cat(" $u\n")
    print(head(object$u[,seq(k)], 2))
    cat("     ...\n\n")

    # cat("Loadings:\n")
    cat(" $v\n")
    print(head(object$v[,seq(k)], 2))
    cat("     ...\n\n")

    # cat("Singular values:\n")
    cat(" $d\n")
    cat(head(object$d, 2))
    cat(" ...\n\n") 
  }
)

#' Print object
#'
#' Print object
#'
#' @param x \code{PCA} object
#' @param ... other arguments
#'
#' @return printed strings
#' @export
#' @rdname print-methods
#' @aliases print,PCA-method
# print.PCA = function(x, ...) {
#   show(x)
# }
setMethod("print", signature(x = "PCA"), 
  function(x, ...) {
  show(x)
})


#' Plot PCAstream
#'
#' Plot PCAstream
#'
#' @param x \code{PCA} object
#' @param main title
#' @param ... other arguments
#'
#' @return plot
#' @export
#' @rdname plot-methods
#' @aliases plot,PCA-method
setMethod("plot", signature(x = "PCA"), 
  function(x, ..., main="scree plot") {

  plot(x$d^2 / x$p, xlab = "Principal component", ylab = "Fraction of total variance", main=main, ...)  
})




#' Evaluate performance of PC estimates
#' 
#' Evaluate performance of PC estimates compared to true PC values
#'
#' @param U true eigen values
#' @param U_est estimated eigen-values
#' @param k number of PCs to evaluate
#' @param metric evaluate the accuracy of the estimated PCs compared to the true PCs using mean explained variance (\code{"MEV"}) or minimum of sum of squared errors (\code{"minSSE"})
#'
#' @return performance metric
#' 
#' @details See performance metrics described by Li, et al. (2023)
#'
#' @references
#' \itemize{
#'  \item Li, Z., Meisner, J., & Albrechtsen, A. (2023). Fast and accurate out-of-core PCA framework for large scale biobank data. Genome Research, 33(9), 1599-1608. \doi{10.1101/gr.277525.122}.
#' }

#' @examples
#' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
#'  X <- hilbert(9)[, 1:6]
#' k <- 4
#' 
#' # Compute SVD using two methods
#' dcmp <- svd( scale(X), k, k)
#' res <- PCAstream( X, k=k)
#' 
#' # Mean variance explained is 1
#' perfMetric(dcmp$u, res$u, metric = "MEV")
#' 
#' # minimum of sum of squared errors zero
#' perfMetric(dcmp$u, res$u, metric = "minSSE")
#
#' @export
perfMetric <- function(U, U_est, k = ncol(U_est), metric = c("MEV", "minSSE")){

  metric <- match.arg( metric )
  stopifnot(is.numeric(k))
  stopifnot(k > 0)
  stopifnot(nrow(U) == nrow(U_est))

  # Modify signs of principal components
  # so diagonal are always positive
  U <- normPC(U[,seq(k),drop=FALSE])
  U_est <- normPC(U_est[,seq(k),drop=FALSE])

  if( metric == "MEV" ){
    values <- vapply(seq(ncol(U_est)), function(j){
      sum(crossprod(U_est, U[,j]))
    }, numeric(1))
    score <- mean(values)
  }else{

    sse <- matrix(NA, ncol(U_est), ncol(U_est))
    for(i in seq(ncol(U_est))){
      for(j in seq(ncol(U_est))){
        sse[i,j] <- sum((U[,j] - U_est[,i])^2)
      }
    }
    score <- sum(apply(sse, 1, min))
  }

  score
}


# Like standard sign function, except sign(x) giving 0 is reset to give 1
sign0 <- function(x) {
  # use standard sign function
  res <- sign(x)

  # get entries that equal 0 and set them to 1
  i <- which(res == 0)
  if (length(i) > 0) {
    res[i] <- 1
  }

  res
}

#' Normalize principal components
#'
#' Modify signs of principal components so diagonal are always positive
#'
#' @param U matrix of principal components
#'
#' @return matrix of principal components with positive diagonal values
#' @examples
#' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
#'  X <- hilbert(9)[, 1:6]
#' k <- 4
#' 
#' dcmp <- svd( scale(X), k, k)
#' 
#' normPC( dcmp$u ) 
#
#' @export
normPC <- function(U){
  sweep(U, 2, sign0(diag(U)), "*")
}









