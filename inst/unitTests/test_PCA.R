test_dashSVD = function(){

	library(RUnit)
	library(GenomicDataStream)

	# simulate n x p matrix with eigen values in d
	simFromSingularValues <- function(n, p, d){
		# 1. Define desired singular values (d) and matrix dimensions (n x p)

		# Ensure the length matches the minimum dimension
		d_vec <- rep(0, min(n, p))
		d_vec[1:length(d)] <- d

		# 2. Create the diagonal matrix D
		D <- diag(d_vec)

		# 3. Create random orthogonal matrices U (n x n) and V (p x p)
		# Using qr.Q(qr(matrix)) ensures orthogonality
		U <- qr.Q(qr(matrix(rnorm(n^2), n, n)))
		V <- qr.Q(qr(matrix(rnorm(p^2), p, p)))

		# 4. Assemble the matrix: A = U * D * t(V)
		A <- U[, 1:min(n, p)] %*% D %*% t(V[, 1:min(n, p)])

		A
	}

	# d.true = sort(c(10, 5, 8, 1), decreasing=TRUE)
	d.true = sort(rexp(10, 1/100), decreasing=TRUE)
	A = simFromSingularValues(1000, 100, d.true)

	k = 10
	dcmp1 <- irlba::irlba(A, k)
	dcmp2 <- GenomicDataStream:::dashSVD_tall(A, k=k)
	dcmp3 <- GenomicDataStream:::dashSVD_wide(A, k=k)

	# check correct singular values
	checkEqualsNumeric( d.true, dcmp1$d[seq(length(d.true))])
	checkEqualsNumeric( d.true, dcmp2$d[seq(length(d.true))])
	checkEqualsNumeric( d.true, dcmp3$d[seq(length(d.true))])

	# check that u is correlated at -1 or 1 with truth
	all(round(diag(crossprod(dcmp1$u, dcmp2$u)), digits=4) %in% c(-1, 1))
	all(round(diag(crossprod(dcmp1$u, dcmp3$u)), digits=4) %in% c(-1,  1))

	# check that v is correlated at -1 or 1 with truth
	all(round(diag(crossprod(dcmp1$v, dcmp2$v)), digits=4) %in% c(-1, 1))
	all(round(diag(crossprod(dcmp1$v, dcmp3$v)), digits=4) %in% c(-1,  1))

}

test_dashSVD_residuals = function(){

	if( require("GenomicDataStreamRegression", quietly=TRUE) ){

		suppressPackageStartupMessages({
		library(DelayedArray)
		library(GenomicDataStreamRegression)
		library(RUnit)
		})

		n <- 2000
		m <- 100
		nc <- 2
		set.seed(1)
		Y <- matrix(rpois(n * m, 100), m, n)
		X <- matrix(rnorm(n * nc), n, nc)
		X <- cbind(1, X) # intercept term

		fit <- glmFitResponses(Y, X, "nb:34")

		A = ResidualMatrixGLM(Y, X, fit = fit, chunkSize = 12)
		Am = as.matrix(A)

		Q = matrix(rnorm(n*m), n,m)

		checkEquals(
		  crossprod(Am, Am %*% Q),
		  crossprod(A, A %*% Q)
		)

		checkEquals(
		  crossprod(Am, Am %*% Q),
		  GenomicDataStream:::mat_prod_AtAQ(Am, Q)
		)

		checkEquals(
		  crossprod(Am, Am %*% Q),
		  GenomicDataStream:::mat_prod_AtAQ(A, Q)
		)
	}
}



test_PCA = function(){

	suppressPackageStartupMessages({
	library(GenomicDataStream)
	library(Matrix)
	library(RUnit)
	library(DelayedArray)})

	# Test on genotype data
	#######################

	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
	files = list.files(dirname(file), "(vcf.gz|bcf|bgen|pgen|bed)$", full.names=TRUE)

	# remove BGEN v1.1, since this doesn't store sample names
	files = grep("1.1.bgen", files, value=TRUE, invert=TRUE)

	k = 4

	for(file in files){
		# cat(file,"\n")

		gds <- GenomicDataStream(file, "GT", init=TRUE) 

		# PCAone
		res <- PCAone(gds, k=k, verbose=FALSE, shuffle=TRUE, threads=1)

		# Standard PCA
		gds <- GenomicDataStream(file, "GT", init=TRUE) 
		dat <- getNextChunk(gds)

		X_scale = scale(dat$X) / sqrt(nrow(dat$X)-1)
		dcmp <- svd(X_scale, k, k)
		
		checkEqualsNumeric( res$d, dcmp$d[1:k], tol=1e-5)
		checkEqualsNumeric( normPC(res$u), normPC(dcmp$u), tol=1e-3)
		checkEqualsNumeric( normPC(res$v), normPC(dcmp$v), tol=1e-3)
	}


	# Test on matrix data
	#####################	

  hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
  X <- hilbert(100)[,1:90]
	k = 4

	# matrix
	X_scale = scale(X) / sqrt(nrow(X) -1)
  dcmp <- irlba::irlba( X_scale, k, k)
  res = dashSVD(X_scale, k)

	checkEqualsNumeric( res$d, dcmp$d[1:k], tol=1e-6)
	checkEqualsNumeric( res$u^2, dcmp$u^2, tol=1e-6)
	checkEqualsNumeric( res$v^2, dcmp$v^2, tol=1e-6)

	# sparseMatrix
	Xa <- as(X, "sparseMatrix")
	X_scale = scale(Xa) / sqrt(nrow(Xa) -1)

  dcmp <- irlba::irlba( X_scale, k, k)
  res = dashSVD(X_scale, k)

	checkEqualsNumeric( res$d, dcmp$d[1:k], tol=1e-6)
	checkEqualsNumeric( res$u^2, dcmp$u^2, tol=1e-6)
	checkEqualsNumeric( res$v^2, dcmp$v^2, tol=1e-6)
	# DelayedArray
	Xa <- DelayedArray(X)

  dcmp <- irlba::irlba( X_scale, k, k)
	# res <- PCAone( Xa, k=k)
	res = dashSVD(DelayedArray(X_scale), k)

	checkEqualsNumeric( res$d, dcmp$d[1:k], tol=1e-6)
	checkEqualsNumeric( res$u^2, dcmp$u^2, tol=1e-6)
	checkEqualsNumeric( res$v^2, dcmp$v^2, tol=1e-6)

	# DelayedArray with transformation
	dcmp <- irlba::irlba( scale(log(X)) / sqrt(nrow(X) -1), k, k)
	# res <- PCAone( log(Xa), k=k)
	res = dashSVD(DelayedArray(scale(log(X)) / sqrt(nrow(X) -1)), k)

	checkEqualsNumeric( res$d, dcmp$d[1:k], tol=1e-6)
	checkEqualsNumeric( res$u^2, dcmp$u^2, tol=1e-6)
	checkEqualsNumeric( res$v^2, dcmp$v^2, tol=1e-6)

	# run SingleCellExperiment
	##########################
	
	library(muscat)
	library(SingleCellExperiment)

	data(example_sce)
	sce <- example_sce

	counts(sce) = DelayedArray(counts(sce))

	# Normalize expression with log2 counts per million
	prior.count <- 1
	lib.size <- colSums2(counts(sce))
	logcounts(sce) <- t(log2(t(counts(sce) + prior.count)) - log2(lib.size) + log2(1e6))

	# PCA on SingleCellExperiment with assay logcounts
	# res <- PCAone(sce, k=k, assay="logcounts")

	X = t(logcounts(sce))
	X_scale = scale(X) / sqrt(nrow(X) -1)
	dcmp = irlba::irlba(X_scale, k, k)

	res = dashSVD(X_scale, k)

	checkEqualsNumeric( res$d, dcmp$d[1:k], tol=1e-3)
	checkEqualsNumeric( res$u^2, dcmp$u^2, tol=1e-3)
	checkEqualsNumeric( res$v^2, dcmp$v^2, tol=1e-3)


	library(beachmat.hdf5)
	library(HDF5Array)

	h5ad_file <- system.file("extdata", "krumsiek11.h5ad",
	                        package="zellkonverter")
	X <- H5ADMatrix(h5ad_file)

	# res <- PCAone(X, k=k)

	X_scale = scale(X) / sqrt(nrow(X) -1)
	dcmp = irlba::irlba(X_scale, k, k)
	res = dashSVD(X_scale, k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( abs(res$u), abs(dcmp$u), tol=1e-3)
	checkEqualsNumeric( abs(res$v), abs(dcmp$v), tol=1e-3)

	if( FALSE ){
		# test on real data
		suppressPackageStartupMessages({
		library(DelayedArray)
		library(SingleCellExperiment)
		library(RUnit)
		library(DelayedMatrixStats)
		library(GenomicDataStream)
		})

		file = "~/Downloads/4e6932db-5a78-40e4-b961-f87f66ba139a.h5ad"

		sce <- readH5AD(file, raw=TRUE)
		sce$total_counts = colSums2(counts(sce))

		lib.size = sce$total_counts
		prior.count = 1
		logcounts(sce) <- t(log2(t(counts(sce) + prior.count)) - log2(lib.size) + log2(1e6))

		k = 10
		system.time({
			res <- dashSVD( logcounts(sce)[seq(100),], k = k)
		})

		system.time({
			dcmp <- irlba::irlba(logcounts(sce)[seq(100),], k, k)
			})

		checkEqualsNumeric(dcmp$d[seq(k)], res$d, tol=1e-3)
		checkEqualsNumeric(abs(dcmp$u), abs(res$u), tol=1e-3)
		checkEqualsNumeric(abs(dcmp$v), abs(res$v), tol=1e-3)
	}

}