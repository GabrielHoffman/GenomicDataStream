
test_DataTable = function(){
	# q()
	# R
	suppressPackageStartupMessages({
	library(GenomicDataStream)
	})

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	file <- system.file("extdata", "test.pvar", package = "GenomicDataStream")
	GenomicDataStream:::test_DataTable( file, "#CHROM")


	file <- system.file("extdata", "test.psam", package = "GenomicDataStream")
	GenomicDataStream:::test_DataTable( file, "#IID")


	file <- system.file("extdata", "test.fam", package = "GenomicDataStream")
	GenomicDataStream:::test_DataTable( file, "")


}





test_pgenstream = function(){

	# q()
	# R
	suppressPackageStartupMessages({
	library(GenomicDataStream)
	})

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	file <- system.file("extdata", "test.pgen", package = "GenomicDataStream")

	# initialize 
	gdsObj = GenomicDataStream(file, chunkSize=3, initialize=TRUE)

	# loop until break
	while( 1 ){

		# get data chunk
		# data$X matrix with features as columns
		# data$info information about each feature as rows
		dat = getNextChunk(gdsObj)

		if( atEndOfStream(gdsObj) ) break
		
		print(dat$info)
	}

     
	file <- system.file("extdata", "test.bed", package = "GenomicDataStream")

	# initialize 
	gdsObj = GenomicDataStream(file, chunkSize=3, initialize=TRUE)

	# loop until break
	while( 1 ){

		# get data chunk
		# data$X matrix with features as columns
		# data$info information about each feature as rows
		dat = getNextChunk(gdsObj)

		if( atEndOfStream(gdsObj) ) break
		
		print(dat$info)
	}


	# Test reading plink 1.x BED
	library(pgenlibr)
	file <- system.file("extdata", "test.bed", package = "GenomicDataStream")

	reg = "1:10001-12000"
	ids = paste0("I", seq(60))
	id_sub = sample(ids, 10)
	ids_str = paste(id_sub, collapse=',')
	idx = sort(match(id_sub, ids))

	gdsObj = GenomicDataStream(file, region=reg, samples=ids_str,chunkSize=3, initialize=TRUE)
	dat = getNextChunk(gdsObj)

	pg = NewPgen(file, raw_sample_ct=60, sample_subset=idx)
	X = ReadList(pg, 2:3)

	max(abs(X - dat$X), na.rm=TRUE)




	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	# GenomicDataStream:::dt(fileIdx)



}




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

	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
	# file <- system.file("extdata", "set1a.vcf.gz", package = "BinaryDosage")
	# system(paste("gunzip -f", file))
	# system(paste("bgzip -f", gsub(".gz$", "", file)))
	# system(paste("tabix -p vcf", file))

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




test = function(){

	# # DEBUG BGEN support on minerva
	# library(GenomicDataStream)
	# library(rbgen)

	# # devtools::reload("~/build2/GenomicDataStream")



	# file = "/hpc/users/hoffmg01/.Rlib/R_433/GenomicDataStream/extdata/test_noph_v1.3_16bits.bgen"
	# samp = c("I1,I2,I3")
	# dat = GenomicDataStream:::getDosage(file, field = '', chunkSize=8, samples=samp)

	# # , samples=c("sample_001","sample_002")
	# rng = data.frame( chromosome = '1', start = 0, end = 12000 )
	# res = bgen.load( file, ranges = rng)

	# samp = c("I1,I2,I3")
	# dat = GenomicDataStream:::getDosage(file, field = '', chunkSize=10, samples=samp)



}


test_large= function(){

	# suppressPackageStartupMessages({
	# library(RUnit)
	# library(VariantAnnotation)
	# library(GenomicDataStream)
	# })

	# # devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	# # Analysis in R
	# #------------------
	# file = "/sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_04032023/data/genotype/PsychAD_genotype.vcf.gz"

	# # VariantAnnotation
	# system.time({
	# gr <- GRanges("chr1", IRanges(0, 1100101))
	# vcf <- suppressWarnings(readVcf(file, param=ScanVcfParam(which=gr)))
	# X_all = geno(vcf)[["DS"]]
	# y = seq(ncol(X_all))

	# res1 = lapply(seq(nrow(X_all)), function(j){
	# 	coef(lm(y ~ X_all[j,]))
	# })
	# res1 = do.call(rbind, res1)
	# rownames(res1) = rownames(X_all)
	# })




	# system.time(res <- GenomicDataStream:::fastLM(y, file, "DS", chunkSize=100000, region="chr1:0-1100001", nthreads=12))








	# X = cbind(1, X_all[1,])
	# summary(lm(y ~ X+0))
	# hatvalues(lm(y ~ X+0))
	# GenomicDataStream:::test_lm(X,y, TRUE)



}


test_xptr = function(){

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	# q()
	# R
	suppressPackageStartupMessages({
	library(GenomicDataStream)
	})

	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
	# file <- system.file("extdata", "test_v1.3_16bits.bgen", package = "GenomicDataStream")


	ptr = GenomicDataStream:::create_xptr(file, "DS", chunkSize=2)

	res = GenomicDataStream:::getNextChunk_rcpp(ptr)


	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

	# initialize 
	gdsObj = GenomicDataStream(file, "DS", chunkSize=5, initialize=TRUE)

	# loop until break
	while( 1 ){

		# get data chunk
		# data$X matrix with features as columns
		# data$info information about each feature as rows
		dat = getNextChunk(gdsObj)

		if( atEndOfStream(gdsObj) ) break
		
		print(dat$info)
	}


		


}




test_chunks = function(){

	# For each file type, create dosage matrix combined across chunks

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

	files = list.files(dirname(file), "(vcf.gz|bcf|bgen)$", full.names=TRUE)

	# test for each file type
	for( file in files){
		X_cat = c()

		# initialize
		obj = GenomicDataStream(file, "DS", chunkSize=10, initialize=TRUE)

		# loop thru chunks
		while( 1 ){

			dat = getNextChunk(obj)

			if( atEndOfStream(obj) ) break
			
			X_cat = cbind(X_cat, dat$X)
		}

		checkEqualsNumeric(t(X_cat), X_all, tol=1e-4) 
	}


}


test_gds_to_fit = function(){

	library(GenomicDataStream)
	library(RUnit)
	library(DelayedArray)

	# lmFitFeatures
	################

	# create response, design and weights 
	y = seq(60)
	X_design = matrix(1, 60,1)
	w = rep(1,60)
	w = seq(60)

	# VCF file
	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

	# Read data into R
	# then run lmFitFeatures()
	gds = GenomicDataStream(file, "DS", initialize=TRUE)
	dat = getNextChunk(gds)
	X_features = dat$X

	res1 = fastlmm::lmFitFeatures(y, X_design, X_features, w)

	# Data stays at C++ level
	# then run lmFitFeatures()
	gds = GenomicDataStream(file, "DS")

	res2 = lmFitFeatures(y, X_design, gds, w)

	sapply(names(res1), function(x){
	 checkEqualsNumeric(res1[[x]], res2[[x]])
	 })

	# R level
	vcf <- suppressWarnings(readVcf(file))
	X_all = geno(vcf)[["DS"]]

	res3 = lapply(seq(nrow(X_all)), function(j){
		coef(lm(y ~ 0 + X_design + X_all[j,], weights=w / mean(w)))
	})
	res3 = do.call(rbind, res3)

	checkEqualsNumeric(res1$coef, res3[,-1], tol=1e-7)


	# lmFitResponses
	################

	n = 100
	m = 55
	nc = 2
	set.seed(1)
	Y = matrix(rnorm(n*m), m, n)
	X = matrix(rnorm(n*nc), n,nc)
	rownames(Y) = seq(m)
	W = matrix(runif(n*m), m,n) 

	fit1 = lmFitResponses(Y, X, W)
	fit2 = lmFitResponses(DelayedArray(Y), X, W)

	sapply(names(fit1), function(x){
	 checkEqualsNumeric(fit1[[x]], fit2[[x]])
	 })

	res3 = lapply(seq(nrow(Y)), function(j){
		coef(lm(Y[j,] ~ 0 + X, weights=W[j,]))
	})
	res3 = do.call(rbind, res3)

	checkEqualsNumeric(fit1$coef, res3, tol=1e-7)


	# trace("lmFitResponses", browser, exit=browser, signature = c("DelayedArray")) 

}

# q()
# R
# library(GenomicDataStream)	
# file = "/Users/gabrielhoffman/prog/R-4.4.1/library/GenomicDataStream/extdata/test.pgen"
# reg = "1:0-150009903"
# gds = GenomicDataStream(file, region=reg, chunkSize=50, initialize=TRUE)
# getNextChunk(gds)$info

# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")
# gds = GenomicDataStream(file, region=reg, chunkSize=2, initialize=TRUE)

# getNextChunk(gds)$info

# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

# ids = "I3,I1,I2"

# reg = "0:0-12000,0:14000-15999"
# gds = GenomicDataStream(file, region=reg, chunkSize=2, samples=ids, initialize=TRUE)

# getNextChunk(gds)


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
	system.time({
	vcf <- suppressWarnings(readVcf(file))
	X_all = geno(vcf)[["DS"]]
	y = seq(ncol(X_all))

	res1 = lapply(seq(nrow(X_all)), function(j){
		coef(lm(y ~ X_all[j,]))
	})
	res1 = do.call(rbind, res1)
	rownames(res1) = rownames(X_all)
	})

	X = cbind(1, X_all[1,])
	# X[,2] = 1
	fit = lm(y ~ X+0)
	fit2 = GenomicDataStream:::test_lm(X,y)
	checkEqualsNumeric(coef(fit), fit2$coefficients)
	
	# Analysis in GenomicDataStream
	#-------------------------------

	files = list.files(dirname(file), "(vcf.gz|bcf|bgen|pgen)$", full.names=TRUE)

	X_design = matrix(1, ncol(X_all))
	w = y
	w[] = 1

	# All regions
	###############
	resList = lapply(files, function(file){

		cat(file, "\n")
		# rm(dat, res)

		# test dosages
		dat = GenomicDataStream:::getDosage(file, "DS", chunkSize=10000)

		# dat$X[1:3, 1:3]
		# t(X_all[1:3,1:3])

		checkEqualsNumeric(t(dat$X), X_all, tol=1e-4) 

		# test regression
		# res <- GenomicDataStream:::lmFitFeatures(y, file, "DS", chunkSize=4)
		gds = GenomicDataStream(file, "DS", chunkSize=5)
		res = lmFitFeatures(y, X_design, gds, w, preprojection=FALSE)

		checkEqualsNumeric(res1, res$coef, silent=TRUE, tol=1e-4)
		})

	# Subset of regions
	####################

	reg = "1:0-150,1:0-15000"

	resList = lapply(files, function(file){

		# cat(file, "\n")

		# test dosages
		dat = GenomicDataStream:::getDosage(file, "DS", region=reg, chunkSize=100)
		checkEqualsNumeric(t(dat$X), X_all[dat$info$ID,], tol=1e-4)

		# test regression
		# res = GenomicDataStream:::lmFitFeatures(y, file, "DS", region=reg)
		gds = GenomicDataStream(file, "DS", region=reg, chunkSize=5)
		res = lmFitFeatures(y, X_design, gds, w, preprojection=FALSE)
		checkEqualsNumeric(res1[dat$info$ID,], res$coef, silent=TRUE, tol=1e-4)
		})


	# Subset of samples
	####################
	reg = "1:0-11000"
	ids = rev(paste0("I", seq(60)))
	ids = sample(ids, 4)
	ids = paste(ids, collapse=',')
	

	# version 1.1 does not support sample ids
	i = grep("v1.1", files)

	resList = lapply(files[-i], function(file){

		cat(file, "\n")

		# test dosages
		dat = GenomicDataStream:::getDosage(file, "DS", region=reg, chunkSize=100, samples=ids)
		checkEqualsNumeric(t(dat$X), X_all[dat$info$ID,rownames(dat$X)], tol=1e-4)
		})

	cbind(dat$X, t(X_all[dat$info$ID,rownames(dat$X),drop=FALSE]))

 	# bcf does not work with subsetting samples
 	# library(vcfppR)
 	# file = "/Users/gabrielhoffman/prog/R-4.4.1/library/GenomicDataStream/extdata/test_noph.bcf"
 	# res = vcftable(file, region="1:0-11000", samples=ids, format="DS")
 	# res$DS

 	# file = "/Users/gabrielhoffman/prog/R-4.4.1/library/GenomicDataStream/extdata/test_noph.vcf.gz"
 	# res = vcftable(file, region="1:0-11000", samples=ids, format="DS")
 	# res$DS


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
		# res = GenomicDataStream:::lmFitFeatures(y, file, "DS", ".")
		gds = GenomicDataStream(file, "DS", region=reg, chunkSize=5)
		res = lmFitFeatures(y, X_design, gds, w, preprojection=FALSE)
		checkEqualsNumeric(res1[rownames(res$coef),], res$coef, silent=TRUE, tol=1e-7)
		})

}


test_DelayedStream = function(){

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")


	suppressPackageStartupMessages({
	library(muscat)
	library(SingleCellExperiment)
	library(Rcpp)
	library(beachmat)
	library(DelayedArray)
	library(RUnit)
	library(GenomicDataStream)
	})

	data(example_sce)

	################
	# sparseMatrix #
	################

	M = counts(example_sce)

	Weights = as.matrix(M)
	Weights[] = 1
	Weights[] = runif(length(Weights))

	design = matrix(1, ncol(M), 1)

	res0 = lapply(seq(nrow(M)), function(j){
		w = Weights[j,]
		w = w / mean(w)
		coef(lm(log(M[j,]+1) ~ 1, weights=w))
	})
	res0 = do.call(rbind, res0)
	rownames(res0) = rownames(M)


	res1 = lmFitResponses(log(M+1), design, Weights)
	checkEqualsNumeric(res1$coef, res0)

	res2 = lmFitResponses(log(as.matrix(M)+1), design, Weights)
	checkEqualsNumeric(res2$coef, res0)

	res3 = lmFitResponses(log(DelayedArray(M)+1), design, Weights)
	checkEqualsNumeric(res3$coef, res0)



	# geneExpr_row = M
	# Weights = as.matrix(geneExpr_row)
	# Weights[] = 1

	# # genes as ROWS
	# res = lmFitResponses(log(geneExpr_row+1), design, Weights)
	# checkEqualsNumeric(res$coef, res1)



	# regrExprResponse = function(Y, chunkSize=114, nthreads=1){

	# 	# wrap wiith beachmat
	# 	ptr = initializeCpp( Y )

	# 	res = GenomicDataStream:::lmFitResponses_export( ptr, design, rownames(Y), chunkSize, nthreads)

	# 	res
	# }

	# res = regrExprResponse(log(M+1))


	# res = GenomicDataStream::lmFitResponses(log(M+1), design, Weights)

	# 	# trace("lmFitResponses", browser, exit=browser, signature = c("dgeMatrix")) 

	
	# checkEqualsNumeric(res$coef, res1)


	# res = GenomicDataStream::lmFitResponses(log(m+1), design, Weights)


	################
	# on-disk h5ad #
	################

	library(RUnit)
	library(SingleCellExperiment)
	library(zellkonverter)

	file <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
	sce <- readH5AD(file, use_hdf5 = TRUE)
	M = assay(sce, 'X')
	Weights = as.matrix(M)
	# Weights[] = 1
	Weights[] = runif(length(Weights))

	res1 = lapply(seq(nrow(M)), function(j){
		w = Weights[j,]
		w = w / mean(w)
		coef(lm(log(M[j,]+1) ~ 1, weights=w))
	})
	res1 = do.call(rbind, res1)
	rownames(res1) = rownames(M)

	design = matrix(1, ncol(M), 1)

	res = lmFitResponses(log(M+1), design, Weights)

	checkEqualsNumeric(res$coef, res1)






}

# test_DelayedStream = function(){

# 	library(GenomicDataStream)

# 	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

# 	# res1 = GenomicDataStream:::getDA_NM( M )

# 	# res2 = GenomicDataStream:::getDA_NM( as.matrix(M) )


# 	library(muscat)
# 	library(SingleCellExperiment)
# 	library(DelayedArray)

# 	data(example_sce)

# 	M = counts(example_sce)[1:3, 1:3]
# 	M = DelayedArray(as.matrix(M))

# 	res1 = GenomicDataStream:::getDA( M )
# 	identical(rownames(M), rownames(res1))
# 	identical(colnames(M), colnames(res1))

# 	res2 = GenomicDataStream:::getDA_NM( as.matrix(M) )
# 	identical(rownames(M), rownames(res2))
# 	identical(colnames(M), colnames(res2))


# 	res1 = GenomicDataStream:::getDA( M )
# 	identical(as.numeric(M), as.numeric(res1))

# 	res2 = GenomicDataStream:::getDA( as.matrix(M) )
# 	identical(as.numeric(M), as.numeric(res2))


# 	res1 = GenomicDataStream:::getDA_eigen( M )
# 	identical(as.numeric(M), as.numeric(res1))

# 	res2 = GenomicDataStream:::getDA_eigen( as.matrix(M) )
# 	identical(as.numeric(M), as.numeric(res2))



# 	res1 = GenomicDataStream:::getDA_vector( M )
# 	identical(as.numeric(M), as.numeric(res1))

# 	res2 = GenomicDataStream:::getDA_vector( as.matrix(M) )
# 	identical(as.numeric(M), as.numeric(res2))


# 	library(Matrix)
# 	library(GenomicDataStream)

# 	x <- round(rsparsematrix(1000, 10, 0.2))

# 	# Initializing it in C++.
# 	library(beachmat)
# 	ptr <- initializeCpp(x)
# 	GenomicDataStream:::column_sums(ptr)
# 	colSums(x)

# }











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



