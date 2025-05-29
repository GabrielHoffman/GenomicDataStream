


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
		res <- PCAstream(gds, k=k, verbose=FALSE, shuffle=FALSE, threads=1)

		# Standard PCA
		gds <- GenomicDataStream(file, "GT", init=TRUE) 
		dat <- getNextChunk(gds)

		X_scale = scale(dat$X) / sqrt(nrow(dat$X)-1)
		dcmp <- svd(X_scale, k, k)
		
		checkEqualsNumeric( res$d, dcmp$d[1:k])
		checkEqualsNumeric( res$u^2, dcmp$u^2)
		checkEqualsNumeric( res$v^2, dcmp$v^2)
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


    hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
    X <- hilbert(9)[, 1:6]
	k = 4

	# matrix
    dcmp <- svd( scale(X), k, k)
	res <- PCAstream( t(X), k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)

	# sparseMatrix
	X <- as(hilbert(9)[, 1:6], "sparseMatrix")

    dcmp <- svd( scale(X), k, k)
	res <- PCAstream( t(X), k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)

	# DelayedArray
	X <- DelayedArray(hilbert(9)[, 1:6])

    dcmp <- svd( scale(X), k, k)
	res <- PCAstream( t(X), k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)

	# DelayedArray with transformation
	dcmp <- svd( scale(log(X)), k, k)
	res <- PCAstream( t(log(X)), k=k)

	checkEqualsNumeric( res$d, dcmp$d[1:k])
	checkEqualsNumeric( res$u^2, dcmp$u^2)
	checkEqualsNumeric( res$v^2, dcmp$v^2)



}