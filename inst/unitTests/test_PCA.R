


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

		# PCAstream
		res <- PCAstream(gds, k=k, verbose=FALSE, shuffle=TRUE, threads=1)

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
  # X <- hilbert(13000)[,seq(10000)]
  X <- hilbert(100)[,1:90]
	k = 4

	# matrix
	X_scale = scale(X) / sqrt(nrow(X) -1)
  dcmp <- irlba::irlba( X_scale, k, k)
	res <- PCAstream( X, k=k, 2, threads=1)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)

	# doesn't give same result due to standardize()
	# res2 <- PCAstream( X, k=k, 100, threads=4, verbose=TRUE)
	# checkEqualsNumeric( res2$d, dcmp$d[1:k])
	# checkEqualsNumeric( res2$u^2, dcmp$u^2)
	# checkEqualsNumeric( res2$v^2, dcmp$v^2)


	# sparseMatrix
	Xa <- as(X, "sparseMatrix")
	X_scale = scale(Xa) / sqrt(nrow(Xa) -1)

  dcmp <- irlba::irlba( X_scale, k, k)
	res <- PCAstream( X, k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)

	# DelayedArray
	Xa <- DelayedArray(X)

  dcmp <- irlba::irlba( X_scale, k, k)
	res <- PCAstream( Xa, k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)

	# DelayedArray with transformation
	dcmp <- irlba::irlba( scale(log(X)) / sqrt(nrow(X) -1), k, k)
	res <- PCAstream( log(Xa), k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)


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
	res <- PCAstream(sce, k=k, assay="logcounts")

	X = t(logcounts(sce))
	X_scale = scale(X) / sqrt(nrow(X) -1)
	dcmp = svd(X_scale, k, k)


	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2, tol=1e-3)
	checkEqualsNumeric( res$v^2, dcmp$v^2, tol=1e-3)











}