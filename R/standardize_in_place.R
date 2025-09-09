


#' Standardize matrix columns in place
#'
#' Standardize mean and variance of matrix columns in place
#' 
#' @param X matrix
#' @param center boolean, TRUE indices center columns
#' @param scale boolean, TRUE indices scale columns
#' 
#' @return none, matrix is standardized in place
#' 
#' @examples
#' X <- matrix(runif(12), 3, 4)
#' standardize_in_place(X)
#' 
#' scale(X)
#' 
#' @export
standardize_in_place <- function(X, center=TRUE, scale=TRUE){

	stopifnot(is.matrix(X))
	stopifnot(is.logical(center))
	stopifnot(is.logical(scale))

	.standardize_in_place(X, center, scale)
}