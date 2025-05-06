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
#'                number of threads (by default \eqn{threads=4}).
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
#' @details PCAstream implements the window-based Randomized SVD proposed by Li et al. 2024

#' @note The singular vectors are not unique and only defined up to sign.
#' If a left singular vector has its sign changed, changing the sign of the corresponding right vector
#' gives an equivalent decomposition.
#'
#' @references
#' \itemize{
#'  \item Z. Li, J Meisner, A Albrechtsen. "Fast and accurate out-of-core PCA framework for large scale biobank data" (2023)
#'        \doi{10.1101/gr.277525.122}.
#' }
#'
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#' 
#' obj <- GenomicDataStream(file, "DS", chunkSize = 3)
#' 
#' res <- PCAstream(obj, k=2)
#' 
#' str(res)
#' 
#' @export
PCAstream <- function(gds, k, p = 7, s = 20, B = 64, threads = 4) {

  stopifnot(is(gds, "GenomicDataStream"))

  chunks <- summaryChunks(gds)
  chunks <- chunks[sample(nrow(chunks)),]  # permute chunks

  N <- slot(gds, "nsamples")
  M <- sum(chunks$counts)

  # M must be large enough, otherwise winSVD is not suggested
  stopifnot(M > B^2)
  # B must be a even
  stopifnot(B %% 2 == 0)
  # nchunks must be larger than B * threads
  stopifnot(nrow(chunks) > B * threads)
  p <- max(c(p, log2(B)+1)) 

  regions <- chunks[sample(nrow(chunks)),"region"]
  region <- paste(regions, collapse = "," )

  res <- stream_pcaone(gds, region, m = M, k = k, s = s, p = p, B = B, threads = threads)
  res
}


#' Get array of chunk sizes
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
  while(1){
    dat <- getNextChunk(x)
    if (atEndOfStream(x)) break
    chunks = rbind(chunks, data.frame(region = paste0(dat$info[1,1], ":", dat$info[1, 2], "-", dat$info[nrow(dat$info),2]), counts = ncol(dat$X)))
  }

  chunks
}




