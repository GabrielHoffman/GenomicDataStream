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
#' 
#' @param B       integer, optional; \cr
#'                number of windows (by default \eqn{B=64}).
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
  H <-  matrix(0, ncol = L,  nrow = N)  ## N x L
  H1 <-  matrix(0, ncol = L,  nrow = N)  ## N x L
  H2 <-  matrix(0, ncol = L,  nrow = N)  ## N x L
  G <- matrix(NA,ncol=L,nrow=M) ## M x L
  switch <- TRUE
  band <- 2
  # block for use
  block <- sample(1:B, M,  replace = T)
  for(i in 1:p) {
    j <- 0
    if (2^(i-1) >= B) {
      H1[] <-  0  ## N x L
      H2[] <-  0  ## N x L
    }
    for(b in 1:B) {
      j <- j + 1
      ## keep <-  (b -(1:band))%%B +1 #used window
      Ab <- A[block==b,]
      Gb <- Ab %*%Omega
      G[block==b,] <- Gb
      if(j <= band / 2) {
        H1 <- H1 + t(Ab)%*%Gb
      } else {
        H2 <- H2 + t(Ab)%*%Gb
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
        H1[] <- 0
        j <- 0
      } else {
        H2[] <- 0
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
  B = qr.solve(t(R), t(H))
  d <- svd(B, nu = k, nv = k)
  list(d = d$d[1:k], u = G %*% d$u, v = d$v)
}

