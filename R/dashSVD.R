
#' Dynamic Shifts Based Randomized SVD
#' 
#' DashSVD for fast randomized SVD
#' 
#' @param x data matrix
#' @param k specifies the target rank of the low-rank decomposition. \eqn{k} should satisfy \eqn{k << min(m,n)}.
#' @param p  number of additional power iterations (by default \eqn{p=7}).
#' @param s oversampling parameter (by default \eqn{s=20}).
#' @param rand randomization method. 1 indicates Gaussian sampling and 2 indicates uniform sampling 
#' @param byrow default \code{FALSE}. Set to true if accessing rows is faster that accessing columns
#' @param verbose print messages
#' 
#' @examples
#' 
#' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
#' X <- hilbert(9)[, 1:6]
#' 
#' dcmp1 <- svd(X)
#' dcmp2 <- dashSVD(X, k=5)
#' 
#' # compare singular values
#' dcmp1$d
#' dcmp2$d
#' 
#' @references
#' Feng, X., Yu, W., Xie, Y. and Tang, J., 2024. Algorithm 1043: Faster randomized svd with dynamic shifts. ACM Transactions on Mathematical Software, 50(2), pp.1-27.
#'
#' @author Zilong Li \email{zilong.dk@gmail.com}, Gabriel Hoffman
#' 
#' @importFrom stats rnorm runif
#' @export
setGeneric(
  "dashSVD",
  function(x, k, p = 3, s = 20, rand = 1, byrow = FALSE, verbose = FALSE) {
    standardGeneric("dashSVD")
})

#' @export
#' @rdname dashSVD
setMethod(
  "dashSVD", signature(x = "ANY"),
  function(x, k, p = 3, s = 20, rand = 1, byrow = FALSE, verbose = FALSE) {

  .dashSVD(x, k = k, p = p, s = s, rand = rand, byrow = byrow, verbose = verbose) 
})

#' @keywords internal
#' @export
.dashSVD <- function(A, k, p = 3, s = 20, rand = 1, byrow = FALSE, verbose = FALSE) {
  if(byrow || nrow(A) > ncol(A)){
    dcmp <- dashSVD_tall(A = A, k = k, p = p, s = s, rand = rand, verbose = verbose)
  } else {
    dcmp <- dashSVD_wide(A = A, k = k, p = p, s = s, rand = rand, verbose = verbose)
  }

  new('PCA', list(
    d = dcmp$d, 
    u = dcmp$u, 
    v = dcmp$v, 
    n = nrow(A), 
    p = ncol(A)))
}

dashSVD_tall <- function(A, k, p = 3, s = 20, rand =  1, verbose = FALSE) {
  M <- nrow(A)
  N <- ncol(A)
  ## L <- k + as.integer(ceiling(k / 2))
  L <- k + s
  if(rand == 1) {
    Omg <- matrix(rnorm(L * M), M, L)
  } else {
    Omg <- matrix(runif(L * M), M, L)
  }
  if( verbose ){
    message("Initial multiplication...")
  }
  Q <- crossprod(A, Omg)
  e <- eigSVD(Q)
  Q <- e$U
  alpha <- 0.0
  for(i in seq(p)) {
    if( verbose){
      message(" power iteration ", i, ", alpha=", round(alpha, 3))
    }
    # e <- eigSVD(t(A) %*% (A %*% Q)-alpha*Q)
    # e <- eigSVD(crossprod(A, A %*% Q) - alpha*Q)
    # mat_prod_AtAQ is faster when A is ResidualMatrixGLM 
    e <- eigSVD(mat_prod_AtAQ(A, Q) - alpha*Q)
    Q <- e$U
    if(e$S[L] > alpha) alpha <- (alpha + e$S[L]) / 2
  }
  if( verbose ){
    message("Final multiplication...")
  }
  e <- eigSVD(A %*% Q)
  U <- e$U
  S <- e$S
  V <- Q %*% e$V

  colnames(U) <- paste0("PC", seq(ncol(U)))
  colnames(V) <- paste0("PC", seq(ncol(V)))

  list(d = S[seq(k)], 
    u = U[,seq(k),drop=FALSE], 
    v = V[,seq(k),drop=FALSE])
}

dashSVD_wide <- function(A, k, p = 3, s = 20, rand =  1, verbose = FALSE) {
  M <- nrow(A)
  N <- ncol(A)
  ## L <- k + as.integer(ceiling(k / 2))
  L <- k + s
  if(rand == 1) {
    Omg <- matrix(rnorm(L * N), N, L)
  } else {
    Omg <- matrix(runif(L * N), N, L)
  }
  if( verbose ){
    message("Initial multiplication...")
  }
  Q <- A %*% Omg
  e <- eigSVD(Q)
  Q <- e$U
  alpha <- 0.0
  for(i in seq(p)) {
    if( verbose){
      message(" power iteration ", i, ", alpha=", round(alpha, 3))
    }
    # e <- eigSVD(   A %*% (t(A) %*% Q) - alpha*Q   )  
    e <- eigSVD(A %*% crossprod(A, Q) - alpha*Q)
    Q <- e$U
    if(e$S[L] > alpha) alpha <- (alpha + e$S[L]) / 2
  }
  if( verbose ){
    message("Final multiplication...")
  }
  e <- eigSVD(crossprod(A, Q))
  V <- e$U
  S <- e$S
  U <- Q %*% e$V

  colnames(U) <- paste0("PC", seq(ncol(U)))
  colnames(V) <- paste0("PC", seq(ncol(V)))

  list(d = S[seq(k)], 
    u = U[,seq(k),drop=FALSE], 
    v = V[,seq(k),drop=FALSE])
}

# SVD via eigendecomposition 
eigSVD <- function(A, tol = 1e-10) {
  m <- nrow(A)
  n <- ncol(A)
  
  # Case 1: Compute via A^T A (smaller covariance matrix)
  eig <- eigen(crossprod(A), symmetric = TRUE)
  V <- eig$vectors
  sigma <- sqrt(pmax(eig$values, 0))  # Ensure non-negative singular values
  
  # Avoid division by near-zero values
  sigma_inv <- ifelse(sigma > tol, 1 / sigma, 0)
  U <- A %*% V %*% diag(sigma_inv)
  
  # Return results as list (thin SVD)
  list(U = U, S = sigma, V = V)
}

# 

#' Evaluate \code{crossprod(A, A \%*\% Q)}
#'
#' When \code{A} is a \code{matrix}, use standard evaluation.  When \code{A} is \code{ResidualMatrixGLM}, use block-wise evaluation 
#'
#' @param A any valid matrix type or \code{ResidualMatrixGLM}
#' @param Q any valid matrix type
#'
#' @rdname mat_prod_AtAQ
#' @keywords internal
#' @export
setGeneric(
  "mat_prod_AtAQ",
  function(A, Q) {
    standardGeneric("mat_prod_AtAQ")
  }
)


#' @rdname mat_prod_AtAQ
#' @keywords internal
#' @export
setMethod(
  "mat_prod_AtAQ", signature(A = "ANY"),
  function(A, Q){
  crossprod(A, A %*% Q)
  }
)

