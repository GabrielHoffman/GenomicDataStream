#' Window-based Randomized SVD
#' 
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
#'                oversampling parameter (by default \eqn{s=10}).
#' 
#' @param B       integer, optional; \cr
#'                number of windows (by default \eqn{B=64}).
#' @param quiet (default: FALSE) suppress messages  
#' 
#' @return \code{winSVD} returns a list containing the following three components:
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
#' @details winSVD implements the window-based Randomized SVD proposed by Li et al. 2024

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
#' @author Zilong Li \email{zilong.dk@gmail.com}, modified by Gabriel Hoffman
#'
#' @examples
#' file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
#' 
#' obj <- GenomicDataStream(file, "DS", chunkSize = 3)
#' 
#' res <- winSVDstream(obj, k=2)
#' 
#' str(res)
#' 
#' @importFrom stats rnorm
#' @importFrom methods slot
#' @importFrom progress progress_bar
#' @export
winSVDstream <- function(gds, k, p = 7, s = 10, B = 64, threads = 16, quiet = FALSE) {

  stopifnot(is(gds, "GenomicDataStream"))  
  chunkSizes <- summaryChunks(gds)
  N <- slot(gds, "nsamples")
  M <- sum(chunkSizes$counts)

  # M must be large enough, otherwise winSVD is not suggested
  stopifnot(M > B^2)
  # B must be a even
  stopifnot(B %% 2 == 0)
  # nchunks must be larger than B * threads
  stopifnot(nrow(chunkSizes) > B * threads)

  p <- max(c(p, log2(B)+1))  
  
  ## each bucket has the index of permuted blocks for each thread
  idx_per_thread <- splitPermuteChunksByThreads(chunkSizes, threads, B)

  if( ! quiet ){
    cat("p: ", p, "\t")
    cat("nchunks: ", nrow(chunkSizes), "\t")
    cat("B: ", B, "\t")
    cat("M: ", M, "\t")
    cat("N: ", N, "\t")
    cat("T: ", length(idx_per_thread), "\n")
  }

  L <- min(k + s, M)
  Omega <- matrix(rnorm(N*L), nrow=N, ncol=L) # N x L
  H <- matrix(0,ncol=L,nrow=N) ## N x L
  H1 <- matrix(0,ncol=L,nrow=N) ## N x L
  H2 <- matrix(0,ncol=L,nrow=N) ## N x L
  G <- matrix(NA, nrow=M, ncol=L) ## M x L

  switch <- TRUE
  band <- 2  ## initial number of buckets as a band

  # indeces were chunks are placed in R
  ## idx = c(1, 1 + cumsum(chunkSizes))

  if (!quiet) {
    pb <- progress_bar$new(
      format = "  svd [:bar] :percent eta: :eta",
      total = p*B, clear = FALSE, width = 60)
  }

  for(i in seq(p)) {
    j <- 0
    if (2^(i-1) >= B) {
      H1[] <-  0  ## N x L
      H2[] <-  0  ## N x L
    }

    ss <- initializeMultiStream(gds, idx_per_thread, chunkSizes)

    bi <- 1
    for(b in seq(B)) {

      if (!quiet) pb$tick()

      j <- j + 1
      nloops <- nrow(idx_per_thread[[1]]) ## the minimum number of loops
      
      for(n in seq(nloops)) {
        Ab <- getChunksParallel(ss[1:threads], threads)
        idx <- unlist(sapply(idx_per_thread[1:threads], function(x) x[n, b]))

        if(length(ss) > threads && n <= nrow(idx_per_thread[[threads+1]]) && !is.na(idx_per_thread[[threads+1]][n,b])) {
          dat <- getNextChunk(ss[[length(ss)]])
          standardize_in_place( dat$X )
          Ab <- cbind( Ab, dat$X )
          idx <- c(idx,  idx_per_thread[[threads+1]][n,b])
        }

        bj <- bi + ncol(Ab) - 1
        if(bj > M) {
          print(Ab)
          print(c(i, bi, bj, ncol(Ab)))
        }
        Gb <- crossprod(Ab, Omega)
        G[seq(bi,bj),] <- Gb
        bi <- bj + 1

        if(j <= band / 2) {
          H1 <- H1 + Ab %*% Gb
        } else {
          H2 <- H2 + Ab %*% Gb
        }
      }

      ## use the first quarter band of succesive iteration (H1)
      ## for extra power iteration updates with the last used band (H2)
      adj <- i>1 & b == 2^(i-2) & 2^(i-1) < B
      if (b < band & !adj ) next
      if (!(j == band || j == band / 2 || adj)) next
      H <- H1 + H2
      QR <- qr(H)
      OmegaOld <- Omega
      Omega <- qr.Q(QR)
      swiched <- colSums(abs(OmegaOld-Omega)) > 2*colSums(abs(OmegaOld+Omega))
      if(switch & any(swiched) & (b+i*2)>3){
        Omega[,swiched] <- -Omega[,swiched]
      }
      if (j == band) {
        H1[] <- 0
        j <- 0
      } else {
        H2[] <- 0
      }
    }
    stopifnot(bj != M)
    band <- min(B,round(band * 2))
  }

  if (!quiet) {
    while (!pb$finished) pb$tick()
    pb$terminate()
  }

  # last step
  getUSV_t(H, G, k)
}

getUSV_t <- function(H,G,k){
  r1 <- qr(G)
  Q <- qr.Q(r1)
  r2 <- qr(Q)
  Rtilt<-  qr.R(r1)
  Rhat <-  qr.R(r2)
  R <- Rhat %*% Rtilt
  B <- qr.solve(t(R), t(H))
  dcmp <- svd(B, nu = k, nv = k)
  list(d = dcmp$d[1:k], u = dcmp$v, v = G %*% dcmp$u)
}

split_into_buckets <- function(sequence, B) {
  n <- length(sequence)
  remainder <- n %% B
  base_size <- floor(n / B)

  if( base_size == 0 ){
    txt = paste0("length of sequence (", n, ") must be larger than B (", B, ")")
    stop(txt)
  }
  
  # Generate bucket indices
  bucket_indices <- rep(seq(B), c(rep(base_size + 1, remainder), rep(base_size, B - remainder)))
  bucket_indices <- sort(bucket_indices)
  
  # Split the sequence
  split(sequence, bucket_indices)
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
  counts <- c()
  chunks <- data.frame(matrix(ncol=3,nrow=0), stringsAsFactors = F)
  # count number of variants
  x <- reinitializeStream(x)
  while(1){
    dat <- getNextChunk(x)
    chunks = rbind(chunks, list(dat$info[1,1], dat$info[1, 2], dat$info[nrow(dat$info),2]+1, ncol(dat$X)))
    if (atEndOfStream(x)) break
  }

  colnames(chunks) <- c("chrom", "start", "end", "counts")
  chunks$region <- paste0(chunks$chrom, ":", chunks$start, "-", chunks$end)
  chunks[,c("region", "counts")]
}

getChunksParallel <- function(gds, ncores) {
  r <- parallel::mclapply(gds, function(x) {
    dat <- getNextChunk(x)
    standardize_in_place( dat$X )
    dat$X
  }, mc.cores = ncores)
  do.call(cbind, r)
}

initializeMultiStream <- function(x, idx_per_thread, chunks) {
  lapply(idx_per_thread, function(idx) {
    idx <- as.vector(idx)
    region <- paste(chunks[idx[!is.na(idx)], "region"], collapse = ",")
    GenomicDataStream(
      file = x@file,
      field = x@field,
      region = region,
      samples = x@samples,
      minVariance = x@minVariance,
      chunkSize = x@chunkSize,
      missingToMean = x@missingToMean,
      initialize = TRUE
    )
  })
}


## regions <- summaryChunks(obj)
## B <- 64
## M <- 6410
## T <- 16

## make a matrix : threads x (nloops x bucket)
splitPermuteChunksByThreads <- function(chunks, threads = 16, B = 64) {
  M <- nrow(chunks)
  idx <- c(rep(1:B, times = M/B), 1:(M %% B))
  chunks_per_window <- split(seq(M), idx)
  T <- threads

  ## w <- chunks_per_window[[1]]
  ## idx <- c(rep(1:T, each = length(w)/T), sample(1:T,length(w) %% T))
  ## chunks_per_thread <- split(w, idx)
  ret <- lapply(1:B, function(b) {
    w <- chunks_per_window[[b]]
    nloops <- length(w)/T
    ## remain <- sample(1:T,length(w) %% T)
    ## idx <- c(rep(1:T, each = nloops), remain)
    idx <- rep(1:T, each = nloops)
    chunks_per_thread <- split(w[1:length(idx)], idx)
    chunks_per_thread[[T+1]] <- w[(length(idx)+1):length(w)]
    chunks_per_thread
  })


  idx_per_thread <- lapply(1:T, function(t) unlist(sapply(ret, function(x) x[[t]])))

  remain <- unlist(sapply(ret, function(x) x[[T+1]] ) )
  if(length(remain)>0) {
    remain <- split_into_buckets(remain, B)
    remain_idx <- t(do.call(rbind, lapply(remain, function(x) {
      length(x) <- max(lengths(remain))  # Set lengths to maximum length
      x
    })))
    idx_per_thread <- c(idx_per_thread, list(remain_idx))
  }

  idx_per_thread

}





