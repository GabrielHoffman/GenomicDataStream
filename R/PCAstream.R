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
winSVDstream <- function(gds, k, p = 7, s = 10, B = 64, quiet = FALSE) {

  stopifnot(is(gds, "GenomicDataStream"))  
    
  gds <- reinitializeStream(gds)
  N <- slot(gds, "nsamples")
  chunkSizes <- countChunks(gds)
  M <- sum(chunkSizes)

  # B must be a even, and larger than length(chunkSizes)
  b <- ceiling(log2(length(chunkSizes)))
  B <- min(B, b + b %% 2)

  stopifnot(B %% 2 == 0)
  p <- max(c(p, log2(B)+1))  
  
  if( ! quiet ){
    cat("B: ", B, "\n")
    cat("p: ", p, "\n")
  }

  L <- min(k + s, M)
  Omega <- matrix(rnorm(N*L), nrow=N, ncol=L) # N x L
  H <-  0  ## N x L
  H1 <-  0  ## N x L
  H2 <-  0  ## N x L
  G <- matrix(NA, nrow=M, ncol=L) ## M x L
  switch <- TRUE

  ## each bucket has the index of blocks
  buckets <- split_into_buckets(seq(length(chunkSizes)), B)
  band <- 2  ## initial number of buckets as a band

  # indeces were chunks are placed in R
  idx = c(1, 1 + cumsum(chunkSizes))

  if (!quiet) {
    pb <- progress_bar$new(
      format = "  svd [:bar] :percent eta: :eta",
      total = p*B, clear = FALSE, width = 60)
  }

  for(i in seq(p)) {
    j <- 0
    if (2^(i-1) >= B) {
      H1 <-  0  ## N x L
      H2 <-  0  ## N x L
    }

    gds <- reinitializeStream(gds)
    
    for(b in seq(B)) {

      if (!quiet) pb$tick()

      j <- j + 1
      for(w in buckets[[b]]) {

        dat <- getNextChunk(gds)
        
        bi <- idx[w]
        bj <- idx[w+1] - 1

        # normalize genotype matrix directly
        # standardize columns in place
        # Ab <- scale(dat$X)
        Ab <- dat$X
        standardize_in_place( Ab )

        Gb <- crossprod(Ab, Omega)
        G[seq(bi,bj),] <- Gb

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
        H1 <- 0
        j <- 0
      } else {
        H2 <- 0
      }
    }
    band <- min(B,round(band * 2))
  }

  if (!quiet) {
    while (!pb$finished) pb$tick()
    pb$terminate()
  }

  # last step
  getUSV_t(H, G, k)
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

