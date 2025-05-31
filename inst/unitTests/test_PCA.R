


test_PCA = function(){

	q()
	R
	suppressPackageStartupMessages({
	library(GenomicDataStream)
	library(Matrix)
	library(RUnit)
	library(DelayedArray)})

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


	# file = files[3]
	# gds <- GenomicDataStream(file, "DS", init=TRUE)
	# a = summaryChunks(gds)
	# gds <- GenomicDataStream(file, "DS", init=TRUE, region=paste0(a$regions, collapse=','))
	# res <- PCAstream(gds, k=k, verbose=TRUE, shuffle=FALSE, p=0)



	# gds <- GenomicDataStream(files[3], "DS", init=TRUE, chunkSize=4)
	# dat <- getNextChunk(gds)
	# str(dat)
	# res <- PCAstream(gds, k=k, verbose=TRUE, shuffle=FALSE, p=0)


	# q()
	# R
	# suppressPackageStartupMessages({
	# library(GenomicDataStream)})

	# file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
	# gds <- GenomicDataStream(file, "DS", init=TRUE)
	# a = summaryChunks(gds)
	# gds <- GenomicDataStream(file, "DS", init=TRUE, region=paste0(a$regions, collapse=','))
	# while (1) {
	#    dat <- getNextChunk(gds)
	 
	#    if (atEndOfStream(gds)) break
	 
	#    print(dat$info)
	#  }


    
	# q()
	# R
	# suppressPackageStartupMessages(
	# library(GenomicDataStream))
	# file = "/Users/gabrielhoffman/prog/R-4.4.2/library/GenomicDataStream/extdata/test.bed"
	# gds1 <- GenomicDataStream(file, "DS", init=TRUE, missingToMean=FALSE)
	# # a = summaryChunks(gds1)
	# dat1 <- getNextChunk(gds1)

	# file = "/Users/gabrielhoffman/prog/R-4.4.2/library/GenomicDataStream/extdata/test.pgen"
	# gds2 <- GenomicDataStream(file, "DS", init=TRUE, missingToMean=FALSE)
	# # summaryChunks(gds2)
	# dat2 <- getNextChunk(gds2)

	# dat1$X[1:3, 1:3]
	# dat2$X[1:3, 1:3]

	# range(dat1$X - dat2$X[,colnames(dat1$X)])

	# i = 6
	# a = dat1$X[, i]
	# b = dat2$X[,colnames(dat1$X)][, i]
	# cbind(a,b)



	q()
	R

	suppressPackageStartupMessages({
	library(GenomicDataStream)
	library(Matrix)
	library(RUnit)
	library(DelayedArray)})

  hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
  # X <- hilbert(13000)[,seq(10000)]
  X <- hilbert(10)
	k = 4

	# matrix
	X_scale = scale(X) / sqrt(nrow(X) -1)
  dcmp <-irlba::irlba( X_scale, k, k)
	res <- PCAstream( t(X), k=k, 2, threads=4, verbose=TRUE)

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
	res <- PCAstream( t(X), k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)

	# DelayedArray
	Xa <- DelayedArray(X)

  dcmp <- irlba::irlba( X_scale, k, k)
	res <- PCAstream( t(Xa), k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)

	# DelayedArray with transformation
	dcmp <- irlba::irlba( scale(log(X)) / sqrt(nrow(X) -1), k, k)
	res <- PCAstream( t(log(Xa)), k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)



}