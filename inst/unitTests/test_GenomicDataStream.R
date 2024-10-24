
test_standardize = function(){

	library(GenomicDataStream)
	library(RUnit)

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	set.seed(1)
	n = 6
	p = 3
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

	suppressPackageStartupMessages({
	library(RUnit)
	library(VariantAnnotation)
	library(vcfppR)
	library(GenomicDataStream)
	})

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	file <- system.file("extdata", "set1a.vcf.gz", package = "BinaryDosage")

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



test_regression = function(){
	
	# Test VCF/BCF/BGEN dosage and regression results

	suppressPackageStartupMessages({
	library(RUnit)
	library(VariantAnnotation)
	library(GenomicDataStream)
	})

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	# Analysis in R
	#------------------
	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

	# VariantAnnotation
	vcf <- suppressWarnings(readVcf(file))
	X_all = geno(vcf)[["DS"]]
	y = seq(ncol(X_all))

	res1 = lapply(seq(nrow(X_all)), function(j){
		coef(lm(y ~ X_all[j,]))
	})
	res1 = do.call(rbind, res1)
	rownames(res1) = rownames(X_all)


	# Analysis in GenomicDataStream
	#-------------------------------

	files = list.files(dirname(file), "(vcf.gz|bcf|bgen)$", full.names=TRUE)

	# All regions
	###############
	resList = lapply(files, function(file){

		# cat(file, "\n")
		# rm(dat, res)

		# test dosages
		dat = GenomicDataStream:::getDosage(file, "DS", chunkSize=100)

		# dat$X[1:3, 1:3]
		# t(X_all[1:3,1:3])

		checkEqualsNumeric(t(dat$X), X_all, tol=1e-4)

		# test regression
		res = GenomicDataStream:::fastLM(y, file, "DS", ".")
		beta = t(do.call(cbind, lapply(res, function(x) x$coef)))
		checkEqualsNumeric(res1, beta, silent=TRUE, tol=1e-4)
		})

	# Subset of regions
	####################

	reg = "1:0-150,1:0-150000"

	resList = lapply(files, function(file){

		# cat(file, "\n")

		# test dosages
		dat = GenomicDataStream:::getDosage(file, "DS", region=reg, chunkSize=100)
		checkEqualsNumeric(t(dat$X), X_all[dat$info$ID,], tol=1e-4)

		# test regression
		res = GenomicDataStream:::fastLM(y, file, "DS", region=reg)
		beta = t(do.call(cbind, lapply(res, function(x) x$coef)))
		checkEqualsNumeric(res1[dat$info$ID,], beta, silent=TRUE, tol=1e-4)
		})


	# Subset of samples
	####################
	reg = "1:0-150,1:0-150000"
	ids = c("I1,I2" )

	# version 1.1 does not support sample ids
	i = grep("v1.1", files)

	resList = lapply(files[-i], function(file){

		# cat(file, "\n")

		# test dosages
		dat = GenomicDataStream:::getDosage(file, "DS", region=reg, chunkSize=100, samples=ids)
		checkEqualsNumeric(t(dat$X), X_all[dat$info$ID,rownames(dat$X)], tol=1e-4)
		})


	# Use GP field from VCF/BCF
	###########################
	files = list.files(dirname(file), "(vcf.gz|bcf)$", full.names=TRUE)

	resList = lapply(files, function(file){

		# cat(file, "\n")
		# rm(dat, res)

		# test dosages
		dat = GenomicDataStream:::getDosage(file, "GP", chunkSize=1000)

		# dat$X[1:3, 1:3]
		# t(X_all[1:3,1:3])

		checkEqualsNumeric(t(dat$X), X_all, tol=1e-4)

		# test regression
		res = GenomicDataStream:::fastLM(y, file, "DS", ".")
		beta = t(do.call(cbind, lapply(res, function(x) x$coef)))
		checkEqualsNumeric(res1, beta, silent=TRUE, tol=1e-4)
		})

}




test_DelayedStream = function(){

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











# test_bgenstream(){



# 	library(RUnit)
# 	library(GenomicDataStream)

# 	file = "/Users/gabrielhoffman/workspace/repos/test_bgen/bgen/example/example.8bits.bgen"

# 	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

# 	# GenomicDataStream:::test_bgen(file)

# 	res = GenomicDataStream:::load( file, sprintf("%s.bgi", file), ranges = data.frame( chromosome = '01', start = 1001, end = 1002 ), rsids = character(0), max_entries_per_sample=3)


# 	head( res$variants )
# 	res$data[1,1:10,1:3]


# 	library( rbgen )

# 	## Test we can load data
# 	D = bgen.load( file, ranges = data.frame( chromosome = '01', start = 0, end = 14444 ), samples=c("sample_001","sample_002"))
# 	# head( D$variants )
# 	# c(D$data)


# 	# D$data[1,1,,drop=FALSE]

# 	# get dosage
# 	X = lapply(seq(nrow(D$data)), function(j){
# 		D$data[j,,] %*% c(0,1,2)
# 		})
# 	names(X) = dimnames(D$data)[[1]]

# 	X = do.call(cbind, X)
# 	colnames(X) = dimnames(D$data)[[1]]


# 	regions = "01:0-14444"
# 	res = GenomicDataStream:::test_bgen( file, "DS", region=regions, samples=c("sample_001,sample_002" ), chunkSize = 1e5)

# 	checkEqualsNumeric(X, res$X)



# 	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream"); rm(res)
# 	regions = "01:0-14444"
# 	res = GenomicDataStream:::test_bgen( file, "DS", region=regions, chunkSize = 1e5, missingToMean=FALSE)

# 	res$X[1:2, 1:2]






# 	regions = "."
# 	y = rnorm(500)


# 	devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream"); rm(res)


# 	res = GenomicDataStream:::test_bgen2( y, file, "DS", region=regions, chunkSize=300)



# 	# Rcpp::sourceCpp(code = "
# 	# 	#include <RcppArmadillo.h>
# 	# 	// [[Rcpp::depends(RcppArmadillo, GenomicDataStream)]]

# 	# 	#include <GenomicDataStream.h>
# 	# 	#include <bgenstream.h>

# 	# 	using namespace std;
# 	# 	using namespace GenomicDataStreamLib;

# 	# 	// [[Rcpp::export]]
# 	# 	string f( string file){
# 	# 		Param param(file, \"DS\");
# 	# 		bgenstream bgenObj(param);
# 	# 		return file;
# 	# 	}")

# 	# f(file)


# 	# Rcpp::cppFunction(code = "
# 	# 	string f( string file){
			
# 	# 	#include <GenomicDataStream.h>
# 	# 		Param param(file, \"DS\");

# 	# 		return file;
# 	# 	}", depends=c("RcppArmadillo", "GenomicDataStream"))

# 	# f(file)

# }



