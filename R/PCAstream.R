#' Window-based Randomized SVD on data streamed from \code{GenomicDataStream}
#' 
#' @param x  \code{GenomicDataStream}
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
#'                number of threads (by default \eqn{threads=4}).  Set to \code{min(threads, floor(nrow(chunks) / B))}
#'
#' @param scaleAndCenter bool, optional; \cr
#'                if \code{TRUE}, scale and center features
#'
#' @param shuffle  bool, optional; \cr
#'                  if \code{TRUE} (default) shuffle genomic regions, the next chunk is not in LD with the previous chunk
#' 
#' @param verbose  string, optional; \cr
#'                  if \code{TRUE} (default is \code{FALSE}) print details
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
#' @details PCAstream implements the window-based Randomized SVD proposed by Li, et al. (2023)

#' @note The singular vectors are not unique and only defined up to sign.
#' If a left singular vector has its sign changed, changing the sign of the corresponding right vector
#' gives an equivalent decomposition.
#'
#' @references
#' \itemize{
#'  \item Li, Z., Meisner, J., & Albrechtsen, A. (2023). Fast and accurate out-of-core PCA framework for large scale biobank data. Genome Research, 33(9), 1599-1608. \doi{10.1101/gr.277525.122}.
#' }
#'
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#' 
#' obj <- GenomicDataStream(file, "DS", chunkSize = 3)
#' 
#' res <- PCAstream(obj, k=5)
#' 
#' res
#' 
#' par(pty="s")
#' plot(res)
#' 
#' @export
setGeneric(
  "PCAstream",
  function(x, k,..., p = 7, s = 20, B = 64, threads = 4, scaleAndCenter = TRUE, shuffle = TRUE, verbose = FALSE) {
    standardGeneric("PCAstream")
  }
)

#' @param algorithm if \code{"serial"} read \code{GenomicDataStream} with a single thread
#'
#' @export
#' @importFrom methods slot
#' @rdname PCAstream
#' @aliases PCAstream,GenomicDataStream-method
setMethod(
  "PCAstream", signature(x = "GenomicDataStream"),
  function(x, k,..., algorithm=c("serial", "parallel"), p = 7, s = 20, B = 64, threads = 4, scaleAndCenter = TRUE, shuffle = TRUE, verbose = FALSE) {

  algorithm <- match.arg(algorithm)
  stopifnot(is(x, "GenomicDataStream"))

  if( verbose ) cat("Read through...\n")

  # get summary of GenomicDataStream
  sObj <- summaryChunks(x)

  # permute chunks
  regions <- sObj$regions

  if( is.null(regions) ){
    stop("Chunks not read from index")
  }

  # if( shuffle ){
  #   regions <- regions[sample(length(regions))]  
  # }

  N <- slot(x, "nsamples")
  M <- sum(sObj$chunks)
  nchunks <- length(sObj$chunks)

  # nchunks must be larger than B * threads
  # stopifnot(nrow(chunks) > B * threads)
  # number of parallel chunks
  n_pll_chunks <- min(threads, floor(length(regions) / B))
  n_pll_chunks <- max(1, n_pll_chunks)

  # set valid B
  while( ! (M > B^2) ){
    B <- B / 2
  }

  # k must be < min(N,M)
  k <- min(c(k,N,M))

  if( verbose ){   
    if( length(N) > 0){ 
      cat(" # samples:", format(N, big.mark=','), "\n")
    }
    cat(" # features:", format(M, big.mark=','), "\n")
    cat(" # chunks:", format(nchunks, big.mark=','), "\n")
    cat(" # threads:", n_pll_chunks, "\n")
    cat(" B:", B, "\n")
    cat(" k:", k, "\n")
  }

  # M must be large enough, otherwise winSVD is not suggested
  # stopifnot(M > B^2)

  # B must be a even
  stopifnot(B %% 2 == 0)

  p <- max(c(p, log2(B)+1)) 

  x <- initializeStream(x)

  if( algorithm == "serial" ){
    # run PCA on GenomicDataStream
    res <- stream_pcaone(x@ptr, 
                      region = paste(regions, collapse = "," ), 
                      m = M, 
                      k = k, 
                      nchunks = nchunks,
                      shuffleRegions = shuffle, 
                      s = s, 
                      p = p, 
                      B = B, 
                      threads = n_pll_chunks, 
                      verbose = verbose)
  }else{
    res <- stream_pcaone_orig(x, 
                      region = paste(regions, collapse = "," ), 
                      m = M, 
                      k = k, 
                      nchunks = nchunks, 
                      s = s, 
                      p = p, 
                      B = B, 
                      threads = n_pll_chunks)
  }

  # set row and column names
  rownames(res$u) <- sObj$sampleIDs
  rownames(res$v) <- sObj$variantIDs

  colnames(res$u) <- paste0("PC", seq(k))
  colnames(res$v) <- paste0("PC", seq(k))
  names(res$d) <- paste0("PC", seq(k))

  new("PCA", res)
})


#' @param chunkSize number of features to read per chunk in \code{GenomicDataStream}
#' 
#' @export
#' @importFrom methods slot
#' @importFrom beachmat initializeCpp
#' @rdname PCAstream
#' @aliases PCAstream,ANY-method
setMethod(
  "PCAstream", signature(x = "ANY"),
  function(x, k, chunkSize = 1000, p = 7, s = 20, B = 64, threads = 4, scaleAndCenter = TRUE, shuffle = TRUE, verbose = FALSE) {

  if( is.null(rownames(x))){
    rownames(x) <- paste0("f_", seq(nrow(x)))
  }

  M <- nrow(x)
  N <- ncol(x)

  # set valid B
  while( ! (M > B^2) ){
    B <- B / 2
  }

  # k must be < min(N,M)
  k <- min(c(k,N,M))

   # B must be a even
  stopifnot(B %% 2 == 0)

  # p <- max(c(p, log2(B)+1)) 

  ptr <- initializeCpp(x)

  res <- stream_pcaone_robj(x = ptr, 
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
                        verbose = verbose)

  # set row and column names
  rownames(res$u) <- colnames(x)
  rownames(res$v) <- rownames(x)

  colnames(res$u) <- paste0("PC", seq(k))
  colnames(res$v) <- paste0("PC", seq(k))
  names(res$d) <- paste0("PC", seq(k))

  new("PCA", res)
})


#' @param chunkSize number of features to read per chunk in \code{GenomicDataStream}
#' @param assay for \code{SummarizedExperiment} which \code{assay} to perform PCA on
#' 
#' @examples
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
  function(x, k, chunkSize = 100, assay="logcounts",..., p = 7, s = 20, B = 64, threads = 4, scaleAndCenter = TRUE, shuffle = TRUE, verbose = FALSE) {

  Y <- SummarizedExperiment::assay(x, assay)
    
  PCAstream(t(Y), k = k, chunkSize = chunkSize, p = p, s = s, B = B, threads = threads, scaleAndCenter = scaleAndCenter, verbose=verbose, ...)
})




#' Get information about chunks
#' 
#' Read through stream to get size of each chunk, IDs of variants and samples
#' 
#' @param x \code{GenomicDataStream}
#' 
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#'
#' # initialize
#' obj <- GenomicDataStream(file, "DS", chunkSize = 3, initialize = TRUE)
#' summaryChunks(obj)
#' 
#' @export
summaryChunks <- function(x){

  x <- reinitializeStream(x)
  summarizeChunks_rcpp(x@ptr)
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
#' @rdname show-methods
#' @importFrom utils head
#' @aliases show,PCA,PCA-method
#' @importFrom methods show
#' @export
setMethod(
  "show", "PCA",
  function(object){ 

    cat("\n       PCA: Computed first", length(object$d), "PCs\n\n")

    k = min(3, length(object$d))

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
#' @param y not used
#' @param ... other arguments
#'
#' @export
#' @rdname plot-methods
#' @aliases plot,PCA-method
setMethod("plot", signature(x = "PCA"), 
  function(x, ...) {

  plot(x$d^2, type="b", ..., xlab = "Principal component", ylab = "Eigen-values")  
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
#' k = 4
#' 
#' dcmp <- svd( scale(X), k, k)
#' res <- PCAstream( t(X), k=k)
#' 
#' perfMetric(dcmp$u, res$u, "MEV")
#' 
#' perfMetric(dcmp$u, res$u, "minSSE")
#
#' @export
perfMetric = function(U, U_est, k = ncol(U_est), metric = c("MEV", "minSSE")){

  metric = match.arg( metric )
  stopifnot(is.numeric(k))

  # Modify signs of principal components
  # so diagonal are always positive
  U <- normPC(U[,seq(k),drop=FALSE])
  U_est <- normPC(U_est[,seq(k),drop=FALSE])

  if( metric == "MEV" ){
    values = sapply(seq(ncol(U_est)), function(j){
      sum(crossprod(U_est, U[,j]))
    })
    score = mean(values)
  }else{

    sse = matrix(NA, ncol(U_est), ncol(U_est))
    for(i in seq(ncol(U_est))){
      for(j in seq(ncol(U_est))){
        sse[i,j] = sum((U[,j] - U_est[,i])^2)
      }
    }
    score = sum(apply(sse, 1, min))
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
#' k = 4
#' 
#' dcmp <- svd( scale(X), k, k)
#' 
#' normPC( dcmp$u ) 
#
#' @export
normPC = function(U){
  sweep(U, 2, sign0(diag(U)), "*")
}









