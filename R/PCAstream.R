#' Window-based Randomized SVD on data streamed from \code{GenomicDataStream}
#' 
#' @param gds  \code{GenomicDataStream}
#'
#' @param k       integer; \cr
#'                specifies the target rank of the low-rank decomposition. \eqn{k} should satisfy \eqn{k << min(m,n)}.
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
#' #' 
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
#' @details PCAstream implements the window-based Randomized SVD proposed by Li, et al. (2024)

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
#' res <- PCAstream(obj, k=5, B=2)
#' 
#' res
#' 
#' par(pty="s")
#' plot(res)
#' 
#' @importFrom methods slot
#' @export
PCAstream <- function(gds, k, p = 7, s = 20, B = 64, threads = 4, shuffle = TRUE, verbose = FALSE) {

  stopifnot(is(gds, "GenomicDataStream"))

  # get summary of GenomicDataStream
  sObj <- summaryChunks(gds)

  # permute chunks
  chunks <- sObj$chunks
  if( shuffle ){
    chunks <- chunks[sample(nrow(chunks)),]  
  }

  N <- slot(gds, "nsamples")
  M <- sum(chunks$counts)

  # nchunks must be larger than B * threads
  # stopifnot(nrow(chunks) > B * threads)
  threads <- min(threads, floor(nrow(chunks) / B))
  threads <- max(1, threads)

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
    cat(" # chunks:", format(nrow(chunks), big.mark=','), "\n")
    cat(" # threads:", threads, "\n")
    cat(" B:", B, "\n")
    cat(" k:", k, "\n")
  }

  # M must be large enough, otherwise winSVD is not suggested
  # stopifnot(M > B^2)

  # B must be a even
  stopifnot(B %% 2 == 0)

  p <- max(c(p, log2(B)+1)) 

  regions <- chunks[sample(nrow(chunks)),"region"]
  region <- paste(regions, collapse = "," )

  gds <- reinitializeStream(gds)

  # run PCA on GenomicDataStream
  res <- stream_pcaone(gds, region, 
                        m = M, 
                        k = k, 
                        s = s, 
                        p = p, 
                        B = B, 
                        threads = threads, 
                        verbose = verbose)

  # set row and column names
  rownames(res$u) <- sObj$sampleIDs
  rownames(res$v) <- sObj$variantIDs

  colnames(res$u) <- paste0("PC", seq(k))
  colnames(res$v) <- paste0("PC", seq(k))
  names(res$d) <- paste0("PC", seq(k))

  # class(res) <- c("PCA", class(res))
  new("PCA", res)
}


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

  chunks <- data.frame(matrix(ncol=2,nrow=0))
  colnames(chunks) <- c("region", "counts")

  # count number of variants
  x <- reinitializeStream(x)
  sampleIDs <- getSampleNames(x)

  lst <- list()
  i <- 1
  while(1){
    dat <- getNextChunk(x)
    if (atEndOfStream(x)) break
    df <- data.frame(region = paste0(dat$info[1,1], ":", dat$info[1, 2], "-", dat$info[nrow(dat$info),2]), counts = ncol(dat$X))
    lst[[i]] <- list(df = df, variantIDs = dat$info$ID)
    i <- i + 1
  }

  chunks <- do.call(rbind, lapply(lst, function(x) x$df))
  variantIDs <- unlist(sapply(lst, function(x) x$variantIDs))

  list( chunks = chunks, 
        sampleIDs = sampleIDs, 
        variantIDs = variantIDs)
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

    cat("Samples:\n")
    cat(" res$u\n")
    print(head(object$u[,seq(k)], 2))
    cat("     ...\n\n")

    cat("Loadings:\n")
    cat(" res$v\n")
    print(head(object$v[,seq(k)], 2))
    cat("     ...\n\n")

    cat("Singular values:\n")
    cat(" res$d\n")
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
#' @aliases print,PCA,PCA-method
# print.PCA = function(x, ...) {
#   show(x)
# }
setMethod("print", signature(x = "PCA"), 
  function(x, ...) {
    cat("PRINT....\n")
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
#' @aliases plot,PCA,PCA-method
setMethod("plot", signature(x = "PCA"), 
  function(x, ...) {

  plot(x$d^2, type="b", ..., xlab = "Principal component", ylab = "Eigen-values")  
})









