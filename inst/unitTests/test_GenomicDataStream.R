
test_DelayedStream(){

	library(GenomicDataStream)

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	# res1 = GenomicDataStream:::getDA_NM( M )

	# res2 = GenomicDataStream:::getDA_NM( as.matrix(M) )


	library(muscat)
	library(SingleCellExperiment)
	library(DelayedArray)

	data(example_sce)

	M = counts(example_sce)[1:3, 1:3]
	M = DelayedArray(as.matrix(M))

	res1 = GenomicDataStream:::getDA_NM( M )
	identical(rownames(M), rownames(res1))
	identical(colnames(M), colnames(res1))

	res2 = GenomicDataStream:::getDA_NM( as.matrix(M) )
	identical(rownames(M), rownames(res2))
	identical(colnames(M), colnames(res2))


	res1 = GenomicDataStream:::getDA( M )
	identical(as.numeric(M), as.numeric(res1))

	res2 = GenomicDataStream:::getDA( as.matrix(M) )
	identical(as.numeric(M), as.numeric(res2))


	res1 = GenomicDataStream:::getDA_eigen( M )
	identical(as.numeric(M), as.numeric(res1))

	res2 = GenomicDataStream:::getDA_eigen( as.matrix(M) )
	identical(as.numeric(M), as.numeric(res2))



	res1 = GenomicDataStream:::getDA_vector( M )
	identical(as.numeric(M), as.numeric(res1))

	res2 = GenomicDataStream:::getDA_vector( as.matrix(M) )
	identical(as.numeric(M), as.numeric(res2))


	library(Matrix)
	library(GenomicDataStream)

	x <- round(rsparsematrix(1000, 10, 0.2))

	# Initializing it in C++.
	library(beachmat)
	ptr <- initializeCpp(x)
	GenomicDataStream:::column_sums(ptr)
	colSums(x)

}


test_standardize = function()

	library(GenomicDataStream)
	library(RUnit)

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	set.seed(1)
	n = 63
	p = 34
	X = matrix(rnorm(n*p), n, p)

	a = base::colSums(X)
	b = GenomicDataStream:::colSums_test(X)
	checkEqualsNumeric(a,b)

	# modify X to make a deep copy

	# TRUE TRUE
	X_res = X + 0.0
	GenomicDataStream:::standardize_test(X_res)
	checkEqualsNumeric(scale(X), X_res)

	# TRUE FALSE
	X_res = X + 0.0
	GenomicDataStream:::standardize_test(X_res, TRUE, FALSE)
	checkEqualsNumeric(scale(X, TRUE, FALSE), X_res)

	# FALSE TRUE
	X_res = X + 0.0
	GenomicDataStream:::standardize_test(X_res, FALSE, TRUE)
	checkEqualsNumeric(scale(X, FALSE, TRUE), X_res)

	# FALSE FALSE
	X_res = X + 0.0
	GenomicDataStream:::standardize_test(X_res, FALSE, FALSE)
	checkEqualsNumeric(scale(X, FALSE, FALSE), X_res)



}



test_vcfstream = function(){

	# q()
	# R

	# file = "/Users/gabrielhoffman/prog/R-4.4.0/library/BinaryDosage/extdata/set2a.vcf.gz"

	# res = GenomicDataStream:::extractVcf( file, "DS", "." )


	suppressPackageStartupMessages({
	library(RUnit)
	library(VariantAnnotation)
	library(vcfppR)
	library(GenomicDataStream)
	})

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	file = "/Users/gabrielhoffman/prog/R-4.4.0/library/BinaryDosage/extdata/set2a.vcf.gz"

	# whole region
	#-------------

	# VariantAnnotation
	vcf <- suppressWarnings(readVcf(file))
	X_all = geno(vcf)[["DS"]]

	res = GenomicDataStream:::extractVcf( file, "DS", "." )
	checkEqualsNumeric(t(res$X), X_all, tol=1e-7)

	res = GenomicDataStream:::extractVcf_eigen( file, "DS", "." )
	checkEqualsNumeric(t(res$X), X_all, tol=1e-7)

	res = GenomicDataStream:::extractVcf_NM( file, "DS", "." )
	checkEqualsNumeric(t(res$X), X_all, tol=1e-7)

	res1 = GenomicDataStream:::extractVcf_vector( file, "DS", "." )
	A = matrix(res1$X, nrow(X_all), ncol(X_all), byrow=TRUE)
	checkEqualsNumeric(as.numeric(A), as.numeric(X_all), tol=1e-7)


	res2 = vcfppR::vcftable(file, ".", format="DS")
	checkEqualsNumeric(res2$DS, X_all, tol=1e-7)

	res3 = GenomicDataStream:::extractVcf_chunks( file, "DS", "." )
	checkEqualsNumeric(res$X[,3:4], res3$X, tol=1e-7)

	# range
	#-------------

	# VariantAnnotation
	region = "1:11000-13000"
	gr = GRanges(1, IRanges(11000, 13000))
	vcf <- suppressWarnings(readVcf(file, "hg19", param=gr))
	X_all = geno(vcf)[["DS"]]

	res = GenomicDataStream:::extractVcf( file, "DS", region )
	checkEqualsNumeric(t(res$X), X_all, tol=1e-7)

	res3 = GenomicDataStream:::extractVcf_chunks( file, "DS", region )
	checkEqualsNumeric(res$X[,3], res3$X, tol=1e-7)

	# subset samples
	#---------------
	idx = c(6, 1, 9, 4)
	sampleIds = paste0(rownames(res$X)[idx], collapse=",")
	res2 = GenomicDataStream:::extractVcf( file, "DS", region, sampleIds )
	checkEqualsNumeric(res2$X, res$X[sort(idx),], tol=1e-7)

	res3 = GenomicDataStream:::extractVcf_chunks( file, "DS", region, sampleIds  )
	checkEqualsNumeric(res$X[sort(idx),3], res3$X, tol=1e-7)


	# DP is an integer
	##################
	file <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
    
    field = "DP"
	vcf <- suppressWarnings(readVcf(file))
	X_all = geno(vcf)[[field]]
	res = GenomicDataStream:::extractVcf( file, field, "." )
	checkEqualsNumeric(t(res$X), X_all)

	res2 = GenomicDataStream:::extractVcf_chunks( file, field, "." )
	checkEqualsNumeric(t(res2$X), X_all[3:4,])

	# GT must be converted to integer
	#################################
	field = "GT"	

	region = "chr21:5030082-5030105"
	res = GenomicDataStream:::extractVcf( file, field, region)
	res2 = GenomicDataStream:::extractVcf_chunks( file, field, region)

	gr = GRanges('chr21', IRanges(5030082, 5030105))
	# vcf <- suppressWarnings(readVcf(file, "hg19", param=gr))
	vcf <- suppressWarnings(readVcf(file, "hg19"))
	X_gt = geno(vcf)[[field]][colnames(res$X),]
	X_ds = X_gt
	X_ds[X_gt=="0/0"] = 0
	X_ds[X_gt=="0/1"] = 1
	X_ds[X_gt=="1/0"] = 1
	X_ds[X_gt=="1/1"] = 2
	X_ds[X_gt=="./."] = NaN
	X_ds = matrix(as.numeric(X_ds), nrow(X_gt), ncol(X_gt), byrow=FALSE)

	checkEqualsNumeric(t(res$X), X_ds[1:3,])
		# X_all[1:3, 1:5]
	# X_all_dosage[1:3, 1:5]
	# t(res$X)[1:3, 1:5]
	checkEqualsNumeric(t(res2$X), X_ds[3,])

	# check missingToMean
	res2 = GenomicDataStream:::extractVcf( file, field, region, missingToMean=TRUE)
	a = colMeans(res$X, na.rm=T)
	b = colMeans(res2$X)
	checkEqualsNumeric(a,b)

	res3 = GenomicDataStream:::extractVcf_chunks( file, field, region, missingToMean=TRUE)
	a = colMeans(res$X, na.rm=T)
	b = colMeans(res3$X)
	checkEqualsNumeric(a[3], b)

	# using multiple chunks
	reg = "chr21:5030082-5030104,chr21:5030105-5030105"
	res4 = GenomicDataStream:::extractVcf_chunks( file, field, reg, missingToMean=TRUE)
	a = colMeans(res$X, na.rm=T)
	b = colMeans(res4$X)
	checkEqualsNumeric(a[3], b)



	# Check errors
	#-------------

	# res = GenomicDataStream:::extractVcf( file, "GT", "." )
	# Error: GT is not supported for multi-allelic site
	# chr21:5030319 chr21:5030319:C:G,T C G,T

	# res = GenomicDataStream:::extractVcf( file, "GT", "343" )
	# Error: region was not found! make sure the region format is correct

     
}


test_regression(){
	
	suppressPackageStartupMessages({
	library(RUnit)
	library(VariantAnnotation)
	library(vcfppR)
	library(GenomicDataStream)
	})

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	file = "/Users/gabrielhoffman/prog/R-4.4.0/library/BinaryDosage/extdata/set2a.vcf.gz"


	# VariantAnnotation
	vcf <- suppressWarnings(readVcf(file))
	X_all = geno(vcf)[["DS"]]
	y = seq(ncol(X_all))

	res1 = lapply(seq(nrow(X_all)), function(j){
		coef(lm(y ~ X_all[j,]))
	})
	res1 = do.call(rbind, res1)

	res = GenomicDataStream:::fastLM(y, file, "DS", "." )
	res = t(do.call(cbind, unlist(res, recursive=FALSE)))
	checkEqualsNumeric(res1, res)

	# add checks removing variants and samples


}





