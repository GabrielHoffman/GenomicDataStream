#' Window-based Randomized SVD
#' 
#' Window-based Randomized SVD on data streamed from \code{GenomicDataStream}
#' 
#' @param A       array_like; \cr
#'                a real/complex \eqn{(m, n)} input matrix (or data frame) to be decomposed.
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
#'                number of windows (by default \eqn{s=64}).
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
#' mat <- matrix(rnorm(100*20000), 20000, 100)
#' res <- winSVD(mat, k = 10)
#' str(res)
#' 
#' @importFrom stats rnorm
#' @importFrom progress progress_bar
#' @export
winSVDstream <- function(gds, k, p = 7, s = 10, B = 64, quiet = FALSE) {
  stopifnot(is(gds, "GenomicDataStream"))

  stopifnot(B %% 2 == 0)
  if(p < log2(B)+1) {
    warning("reset p as log2(B)+1, which is ", log2(B)+1)
  }
  p <- max(c(p, log2(B)+1))
  
  # count number of variants
  M <- 0
  gds <- reinitializeStream(gds)
  while(1){
    dat <- getNextChunk(gds)
    if (atEndOfStream(gds)) break
  }
  M = featuresRead(gds)
  message("# variants:", M)
    
  N <- slot(gds, "nsamples")
  cs <- slot(gds, "chunkSize")
  
  L <- k + s
  Omega <- matrix(rnorm(N*L),nrow=N,ncol=L) # N x L
  H <- matrix(0,ncol=L,nrow=N) ## N x L
  H1 <- matrix(0,ncol=L,nrow=N) ## N x L
  H2 <- matrix(0,ncol=L,nrow=N) ## N x L
  G <- matrix(NA,ncol=L,nrow=M) ## M x L
  switch <- TRUE
  ## how many times we call getNextChunk
  nchunks <- ceiling(M / cs)
  ## each bucket has the index of blocks
  buckets <- split_into_buckets(seq(nchunks), B)
  band <- 2  ## initial number of buckets as a band

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

    gds <- reinitializeStream(gds)
    
    for(b in seq(B)) {

      if (!quiet) pb$tick()

      j <- j + 1
      bi <- 1
      for(w in buckets[[b]]) {
        dat <- getNextChunk(gds)
        bi <- (w-1) * ncol(dat$X) + 1
        ## normalize genotype matrix directly
        Ab <- scale(dat$X)

        bj <- bi + ncol(Ab) - 1
        if(bj > M) bj <- M

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
      adj <-  i>1 & b == 2^(i-2) & 2^(i-1) < B
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
    band <- min(B,round(band * 2))
  }

  if (!quiet) {
    while (!pb$finished) pb$tick()
    pb$terminate()
  }

  # keep rows with valid entries
  # tail can have NA values since each chunk has a cap of chunkSize
  # but likely has fewer variants after MAF filter
  keep <- !is.na(G[,1])

  # last step
  getUSV(H, G[keep,], k)
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

# a = apply(dat$X, 2, maf)
# keep = (a > MAF)
# Ab <- scale(dat$X)[,keep]

maf = function(x){
  p = sum(x) / (2*length(x))
  min(p, 1-p)
}



