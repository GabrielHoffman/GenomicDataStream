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
#'                number of additional power iterations (by default \eqn{p=8}).
#'
#' @param s       integer, optional; \cr
#'                oversampling parameter (by default \eqn{s=20}).
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
winSVDstream <- function(gds, k, p = 8, s = 20, B = 64, threads = 4, quiet = FALSE) {

  stopifnot(is(gds, "GenomicDataStream"))  
  chunks <- summaryChunks(gds)
  N <- slot(gds, "nsamples")
  M <- sum(chunks$counts)

  # M must be large enough, otherwise winSVD is not suggested
  stopifnot(M > B^2)
  # B must be a even
  stopifnot(B %% 2 == 0)
  # nchunks must be larger than B * threads
  stopifnot(nrow(chunks) > B * threads)

  p <- max(c(p, log2(B)+1))  
  
  ## each bucket has the index of permuted blocks for each thread
  idx_per_thread <- splitPermuteChunksByThreads(chunks, threads, B)

  if( ! quiet ){
    cat("p: ", p, "\t")
    cat("nchunks: ", nrow(chunks), "\t")
    cat("B: ", B, "\t")
    cat("M: ", M, "\t")
    cat("N: ", N, "\t")
    cat("T: ", length(idx_per_thread), "\n")
  }

  L <- min(k + s, M/2)
  Omega <- matrix(rnorm(N*L), nrow=N, ncol=L) # N x L
  H <- matrix(0,ncol=L,nrow=N) ## N x L
  H1 <- matrix(0,ncol=L,nrow=N) ## N x L
  H2 <- matrix(0,ncol=L,nrow=N) ## N x L
  G <- matrix(NA, nrow=M, ncol=L) ## M x L

  switch <- TRUE
  band <- 2  ## initial number of buckets as a band

  if (!quiet) {
    pb <- progress_bar$new(
      format = "  svd [:bar] :percent eta: :eta",
      total = p*B, clear = FALSE, width = 60)
  }

  ids <- c() ## SNP ids
  
  for(i in seq(p)) {
    j <- 0
    if (2^(i-1) >= B) {
      H1[] <-  0  ## N x L
      H2[] <-  0  ## N x L
    }


    nloops <- sum(complete.cases(idx_per_thread[[1]])) ## the minimum number of loops
    extra_thread <- NULL
    if(length(idx_per_thread) > threads || threads == 1) extra_thread <- idx_per_thread[[length(idx_per_thread)]]
    bi <- 1

    for(b in seq(B)) {

      if (!quiet) pb$tick()

      j <- j + 1
      
      for(n in seq(nloops)) {
        #TODO: permute the order of columns of Ab
        Ab <- getChunksParallel(gds, idx_per_thread[seq(threads)], chunks, n, b)

        if(!is.null(extra_thread ) && n <= nrow(extra_thread) && ncol(extra_thread) >= b && !is.na(extra_thread[n,b])) {
          regions <- chunks[as.vector(na.omit(extra_thread[n:nrow(extra_thread),b])),"region"]
          regions <- paste(regions, collapse = ",")
          ss <- reinitializeStream(gds, region = regions)
          if(threads > 1) {
            dat <- getNextChunk(ss)
            standardize_in_place( dat$X )
            Ab <- cbind( Ab, dat$X )
          }

          # in case nloops can be smaller than nrow of extra_thread
          if(n == nloops && nloops < nrow(extra_thread)) {
            for(n in seq(nloops+1, nrow(extra_thread))) {
              if(is.na(extra_thread[n,b])) break
              dat <- getNextChunk(ss)
              standardize_in_place( dat$X )
              Ab <- cbind( Ab, dat$X )
            }
          }
        }

        if(i == p) ids <- c(ids, colnames(Ab)) ## store ids in the last epoch
        
        bj <- bi + ncol(Ab) - 1
        ## print(c(bi, bj, ncol(Ab), n, b))
        stopifnot(bj < M + 1)
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

    stopifnot(bj == M)
    band <- min(B,round(band * 2))
  }

  if (!quiet) {
    while (!pb$finished) pb$tick()
    pb$terminate()
  }

  # last step
  ret <- getUSV_t(H, G, k)
  rownames(ret$v) <- ids
  return(ret)
}

getUSV_t <- function(H, G, k){
  r1 <- qr(G)
  Q <- qr.Q(r1)
  r2 <- qr(Q)
  Rtilt<-  qr.R(r1)
  Rhat <-  qr.R(r2)
  R <- Rhat %*% Rtilt
  B = qr.solve(t(R), t(H))
  d <- svd(B, nu = k, nv = k)

  list(d = d$d[1:k], u = d$v, v = G %*% d$u)
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

getChunksParallel <- function(x, idx_per_thread, chunks, n, b) {
  res <- parallel::mclapply(idx_per_thread, function(idx) {
    if(is.na(idx[n,b])) return(NULL)
    region <- chunks[idx[n, b], "region"]
    gds <- GenomicDataStream(
      file = x@file,
      field = x@field,
      region = region,
      samples = x@samples,
      minVariance = x@minVariance,
      chunkSize = x@chunkSize,
      missingToMean = x@missingToMean,
      initialize = TRUE
    )
    if (atEndOfStream(gds)) return(NULL)
    dat <- getNextChunk(gds)
    standardize_in_place( dat$X )
    return(dat$X)
  }, mc.cores = length(idx_per_thread))

  do.call(cbind, res)
}


## make a matrix : threads x (nloops x bucket)
splitPermuteChunksByThreads <- function(chunks, threads = 16, B = 64) {

  M <- nrow(chunks)

  idx <- rep(1:B, times = ceiling(M/B))
  chunks_per_window <- split(c(seq(M), rep(NA, length(idx)-M)), idx)

  ret <- lapply(1:B, function(b) {
    w <- chunks_per_window[[b]]
    nloops <- floor(length(w)/threads)
    idx <- rep(1:threads, each = nloops)
    chunks_per_thread <- split(w[1:length(idx)], idx)
    chunks_per_thread
  })


  idx_per_thread <- lapply(seq(threads), function(t) unlist(sapply(ret, function(x) x[[t]])))

  remain <- setdiff(seq(M), as.vector(unlist(ret)))

  if(length(remain) > 0) {
    idx <- rep(1:B, times = ceiling(length(remain)/B))

    remain_idx <- split(c(remain, rep(NA, length(idx)-length(remain))), idx)

    remain_idx <- (do.call(cbind, lapply(remain_idx, function(x) {
      length(x) <- max(lengths(remain_idx))  # Set lengths to maximum length
      x
    })))

    if(ncol(remain_idx) == B && nrow(remain_idx) > 1) {
      # the last chunk may have less variants than specified chunkSize
      # we want to read it at the final window, otherwise the actual chunksize may not match input chunkSize
      i <- tail(which(!is.na(remain_idx[,B])), 1) ## last non na index

      last <- remain_idx[i,B]
      remain_idx[which(remain_idx == M)] <- last
      remain_idx[i, B] <- M
    }

    idx_per_thread <- c(idx_per_thread, list(remain_idx))
  }

  if(threads == 1) {
    a <- idx_per_thread[[1]]
    i <- tail(which(!is.na(a[,B])), 1) ## last non na index
    last <- a[i,B]
    a[which(a == M)] <- last
    a[i, B] <- M
    idx_per_thread[[1]] <- a
  }

  idx_per_thread

}





