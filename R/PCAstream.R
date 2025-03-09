#' @note assume the input matrix  is tall, i.e. nrow (features) > ncol (samples)
#' @export
winSVDstream <- function(gds, meta, k, p = 7, s = 10, B = 64) {
  stopifnot(is(gds, "GenomicDataStream"))
  stopifnot(is(meta, "metadata"))
  stopifnot(B %% 2 == 0)
  if(p < log2(B)+1) {
    warning("reset p as log2(B)+1, which is ", log2(B)+1)
  }
  p <- max(c(p, log2(B)+1))
  
  N <- meta$N
  M <- meta$M
  cs <- gds$chunkSize
  
  L <- k + s
  Omega <- matrix(stats::rnorm(N*L),nrow=N,ncol=L) # N x L
  H <-  0  ## N x L
  H1 <-  0  ## N x L
  H2 <-  0  ## N x L
  G <- matrix(NA,ncol=L,nrow=M) ## M x L
  switch <- TRUE
  nchunks <- ceiling(M / cs) ## how many times we call getNextChunk
  ## each bucket has the index of blocks
  buckets <- split_into_buckets(1:nchunks, B)
  band <- 2  ## initial number of buckets as a band
  
  for(i in 1:p) {
    j <- 0
    if (2^(i-1) >= B) {
      H1 <-  0  ## N x L
      H2 <-  0  ## N x L
    }
    
    for(b in 1:B) {
      j <- j + 1
      for(w in buckets[[b]]) {
        dat <- getNextChunk(obj)
        bi <- (w-1) * cs + 1
        bj <- w * cs
        if(bj > M) bj <- M
        Ab <- (dat$X - meta$mean[bi:bj]) * meta$scale[bi:bj]  ## normalize
        Gb <- Ab %*%Omega
        G[bi:bj,] <- Gb
        if(j <= band / 2) {
          H1 <- H1 + t(Ab)%*%Gb
        } else {
          H2 <- H2 + t(Ab)%*%Gb
        }
      }
      ## use the first quarter band of succesive iteration (H1)
      ## for extra power iteration updates with the last used band (H2)
      adj <-  i>1 & b == 2^(i-2) & 2^(i-1) < B
      if (b < band & !adj ) next
      if (!(j == band || j == band / 2 || adj)) next
      H <- H1 + H2
      QR<-qr(H)
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
  getUSV(H, G, k)
}

## for stream version we need to know
## 1. the dimension of data
## 2. how to normalize the data
#' @export
PCAstream_prepare <- function(gds, center = TRUE, scale = TRUE) {
  stopifnot(is(gds, "GenomicDataStream"))
  res <- list(M = 0, N = 0, mean = NULL, scale = NULL)
  # loop until break
  while (1) {
    # get data chunk
    # data$X matrix with features as columns
    # data$info information about each feature as rows
    dat <- getNextChunk(gds)
    if (atEndOfStream(gds)) break
    if(res$M > 0) stopifnot(res$N == nrow(dat$X))
    res$M <- res$M + ncol(dat$X)
    res$N <- nrow(dat$X)
    res$mean <- c(res$mean, colMeans(dat$X))
    scale <- apply(dat$X, 2, sd)
    scale[scale <= 0] <- 1 ## in case of any fixed SNP
    res$scale <- c(res$scale, 1.0 / scale) ## multiply this is faster
  }
  class(res) <- "metadata"
  res
}

split_into_buckets <- function(sequence, B) {
  n <- length(sequence)
  remainder <- n %% B
  base_size <- floor(n / B)
  
  # Generate bucket indices
  bucket_indices <- rep(1:B, c(rep(base_size + 1, remainder), rep(base_size, B - remainder)))
  bucket_indices <- sort(bucket_indices)
  
  # Split the sequence
  split(sequence, bucket_indices)
}

