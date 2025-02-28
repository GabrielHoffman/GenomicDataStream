#' @title window-based Randomized SVD
#' 
#' @description winSVD implements the window-based Randomized SVD proposed by Li et al. 2024
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
#' @note The singular vectors are not unique and only defined up to sign.
#' If a left singular vector has its sign changed, changing the sign of the corresponding right vector
#' gives an equivalent decomposition.
#'
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
#' mat <- matrix(rnorm(100*20000), 20000, 100)
#' res <- winSVD(mat, k = 10)
#' str(res)
#' @export
winSVD <- function(A, k, p = 7, s = 10, B = 64) {
  stopifnot(B %% 2 == 0)
  if(p < log2(B)+1) {
    warning("reset p as log2(B)+1, which is ", log2(B)+1)
  }
  p <- max(c(p, log2(B)+1))
  if(nrow(A) < ncol(A)) A <- t(A) ## make it tall 
  N <- ncol(A)
  M <- nrow(A)
  L <- k + s
  Omega <- matrix(stats::rnorm(N*L),nrow=N,ncol=L) # N x L
  H <-  0  ## N x L
  H1 <-  0  ## N x L
  H2 <-  0  ## N x L
  G <- matrix(NA,ncol=L,nrow=M) ## M x L
  switch <- TRUE
  band <- 2
  # block for use
  block <- sample(1:B, M,  replace = T)
  for(i in 1:p) {
    j <- 0
    if (2^i >= B) {
      H1 <-  0  ## N x L
      H2 <-  0  ## N x L
    }
    for(b in 1:B) {
      j <- j + 1
      ## keep <-  (b -(1:band))%%B +1 #used window
      Ab <- A[block==b,]
      Gb <- Ab %*%Omega
      G[block==b,] <- Gb
      if(j <= band / 2) {
        H1 <- H1 + t(Ab)%*%Gb
      } else if (j <= band) {
        H2 <- H2 + t(Ab)%*%Gb
      }

      if((b-1)%%2^(i-1) != 2^(i-1)-1 & b!=B)
        next
      
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
      }
      if (j == band / 2) {
        H2 <- 0
      }
    }
    band <- min(B,round(band * 2))

  }
  getUSV(H, G, k)
}

getUSV <- function(H,G,k){
  r1 <- qr(G)
  Q <- qr.Q(r1)
  r2 <- qr(Q)
  Rtilt<-  qr.R(r1)
  Rhat <-  qr.R(r2)
  R <- Rhat%*%Rtilt
  ## R.T * B = H.T
  ## B <- t(MASS::ginv(R)) %*% t(H)
  B = qr.solve(t(R), t(H))
  d <- svd(B, nu = k, nv = k)
  list(d = d$d[1:k], u = G %*% d$u, v = d$v)
}

## for stream version we need to know
## 1. the dimension of data
## 2. how to normalize the data
streamPCA_prepare <- function(gds, center = TRUE, scale = TRUE) {
  res <- list(M = 0, N = 0, C = NULL, S = NULL)
  # loop until break
  while (1) {
    # get data chunk
    # data$X matrix with features as columns
    # data$info information about each feature as rows
    dat <- getNextChunk(obj)
    if(res$M > 0) stopifnot(res$N == nrow(dat$X))
    res$M <- res$M + ncol(dat$X)
    res$N <- nrow(dat$X)
    avg <- colMeans(dat$X)
    res$C <- c(res$C, avg)
    res$S <- c(res$S, 1.0/(avg*(1.0-avg))) # for genetic data
    if (atEndOfStream(obj)) break
  }
  res
}
