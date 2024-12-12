

test = function(){

	library(RUnit)

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

	res = GenomicDataStream:::f_test_mat_types( file )

	sapply(res, function(x){
		checkEqualsNumeric(x, res[[1]])
		})




}




test_DataTable = function(){
	# q()
	# R
	suppressPackageStartupMessages({
	library(GenomicDataStream)
	})

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	file <- system.file("extdata", "test.pvar", package = "GenomicDataStream")
	res1 = GenomicDataStream:::test_DataTable( file, "#CHROM")


	file <- system.file("extdata", "test.psam", package = "GenomicDataStream")
	res2 = GenomicDataStream:::test_DataTable( file, "#IID")


	file <- system.file("extdata", "test.fam", package = "GenomicDataStream")
	res3 = GenomicDataStream:::test_DataTable( file, "")


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
		
		# print(dat$info)
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
		
		# print(dat$info)
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
		gds = GenomicDataStream(file, "DS", chunkSize=10000, initialize=TRUE, missingToMean=TRUE)
	    dat = getNextChunk(gds)

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
		gds = GenomicDataStream(file, "DS", region=reg, chunkSize=10000, initialize=TRUE, missingToMean=TRUE)
	    dat = getNextChunk(gds)
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
	ids = sample(ids, 40)
	ids = paste(ids, collapse=',')
	

	# version 1.1 does not support sample ids
	i = grep("v1.1", files)

	resList = lapply(files[-i], function(file){

		cat(file, "\n")

		# test dosages
		# dat = GenomicDataStream:::getDosage(file, "DS", region=reg, chunkSize=100, samples=ids)

		gds = GenomicDataStream(file, "DS", region=reg, chunkSize=10000, initialize=TRUE, samples=ids, missingToMean=TRUE)
	    dat = getNextChunk(gds)
		checkEqualsNumeric(t(dat$X), X_all[dat$info$ID,rownames(dat$X)], tol=1e-3)
		})

	# cbind(dat$X, t(X_all[dat$info$ID,rownames(dat$X),drop=FALSE]))

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
		gds = GenomicDataStream(file, "GP", chunkSize=10000, initialize=TRUE)
	    dat = getNextChunk(gds)

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


test_glmFitFeatures = function(){


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

	X_design = matrix(1, ncol(X_all))
	w = y
	w[] = 1
	w = seq(length(w))
	offset = log(seq(length(y)))

	# GAUSSIAN
	##########
	set.seed(1)
	y = rnorm(nrow(X_design))
	gds = GenomicDataStream(file, "DS")
	res1 = glmFitFeatures(y, X_design, gds, "gaussian", w, offset)

	gds = GenomicDataStream(file, "DS", initialize=TRUE)
	dat <- getNextChunk(gds)
	res2 = glmFitFeatures(y, X_design, dat$X, "gaussian", w, offset)

	sapply(names(res1), function(x){
		checkEqualsNumeric(res1[[x]], res2[[x]])
		})

	# LOGIT
	#######
	set.seed(1)
	y = sample(2, nrow(X_design), replace=TRUE) - 1
	gds = GenomicDataStream(file, "DS")
	res1 = glmFitFeatures(y, X_design, gds, "logit", w, offset)

	gds = GenomicDataStream(file, "DS", initialize=TRUE)
	dat <- getNextChunk(gds)
	res2 = glmFitFeatures(y, X_design, dat$X, "logit", w, offset)

	sapply(names(res1), function(x){
		checkEqualsNumeric(res1[[x]], res2[[x]])
		})

	# PROBIT
	########

	gds = GenomicDataStream(file, "DS")
	res1 = glmFitFeatures(y, X_design, gds, "probit", w, offset)

	gds = GenomicDataStream(file, "DS", initialize=TRUE)
	dat <- getNextChunk(gds)
	res2 = glmFitFeatures(y, X_design, dat$X, "probit", w, offset)

	sapply(names(res1), function(x){
		checkEqualsNumeric(res1[[x]], res2[[x]])
		})

	# POISSON
	##########

	set.seed(1)
	y = rpois(nrow(X_design), 100)
	gds = GenomicDataStream(file, "DS")
	res1 = glmFitFeatures(y, X_design, gds, "poisson", w, offset)

	gds = GenomicDataStream(file, "DS", initialize=TRUE)
	dat <- getNextChunk(gds)
	res2 = glmFitFeatures(y, X_design, dat$X, "poisson", w, offset)

	sapply(names(res1), function(x){
		checkEqualsNumeric(res1[[x]], res2[[x]])
		})

	# NB:10
	#######

	gds = GenomicDataStream(file, "DS")
	res1 = glmFitFeatures(y, X_design, gds, "nb:10", w, offset)

	gds = GenomicDataStream(file, "DS", initialize=TRUE)
	dat <- getNextChunk(gds)
	res2 = glmFitFeatures(y, X_design, dat$X, "nb:10", w, offset)

	sapply(names(res1), function(x){
		checkEqualsNumeric(res1[[x]], res2[[x]])
		})

	# NB
	#####

	gds = GenomicDataStream(file, "DS")
	res1 = glmFitFeatures(y, X_design, gds, "nb", w, offset)

	gds = GenomicDataStream(file, "DS", initialize=TRUE)
	dat <- getNextChunk(gds)
	res2 = glmFitFeatures(y, X_design, dat$X, "nb", w, offset)

	sapply(names(res1), function(x){
		checkEqualsNumeric(res1[[x]], res2[[x]])
		})




}

test_glmFitResponses = function(){

	suppressPackageStartupMessages({
	library(RUnit)
	library(GenomicDataStream)
	library(DelayedArray)
	library(MASS)
	})

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	n = 1000
	m = 100
	nc = 2
	set.seed(1)
	Y = matrix(rnorm(n*m), m, n)
	colnames(Y) = paste0("s", seq(n))
	rownames(Y) = paste0("r", seq(m))

	X = cbind(1, matrix(rnorm(n*nc), n,nc))
	colnames(X) = paste0("V", seq(nc+1))
	W = matrix(runif(n*m), ,n)	

	isSame = function(res1, res2){
		sapply(names(res1), function(x){
		checkEqualsNumeric(res1[[x]], res2[[x]])
		})	
	}

	# GAUSSIAN
	res0 = lmFitResponses(Y, X, detail=4)
	res1 = fastLinReg::glmFitResponses(Y, X, "gaussian", detail=4)
	res2 = GenomicDataStream::glmFitResponses(Y, X, "gaussian", detail=4)
	isSame(res1, res2)

	# PROBIT
	for(i in seq(nrow(Y))){
		Y[i,] = sample(2, n, replace=TRUE) - 1
	}

	res1 = fastLinReg::glmFitResponses(Y, X, "probit")
	res2 = GenomicDataStream::glmFitResponses(Y, X, "probit")
	isSame(res1, res2)

	# LOGIT
	res1 = fastLinReg::glmFitResponses(Y, X, "logit")
	res2 = GenomicDataStream::glmFitResponses(Y, X, "logit")
	isSame(res1, res2)

	# POISSON
	for(i in seq(nrow(Y))){
		Y[i,] = rnegbin(n, 5, 3) 
	}
	
	res1 = fastLinReg::glmFitResponses(Y, X, "poisson")
	res2 = GenomicDataStream::glmFitResponses(Y, X, "poisson")
	isSame(res1, res2)

	# NB
	res1 = fastLinReg::glmFitResponses(Y, X, "nb")
	res2 = GenomicDataStream::glmFitResponses(Y, X, "nb")
	isSame(res1, res2)

	res1 = GenomicDataStream::glmFitResponses(DelayedArray(Y), X, "gaussian")
	res2 = GenomicDataStream::glmFitResponses(DelayedArray(Y), X, "gaussian")
	isSame(res1, res2)


	# q()
	# R
	# suppressPackageStartupMessages({
	# library(RUnit)
	# library(GenomicDataStream)
	# library(DelayedArray)
	# library(MASS)
	# }) 

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	n = 1000
	m = 100
	nc = 2
	set.seed(1)
	Y = matrix(rnorm(n*m), m, n)
	colnames(Y) = paste0("s", seq(n))
	rownames(Y) = paste0("r", seq(m))
	X = cbind(1, matrix(rnorm(n*nc), n,nc))
	colnames(X) = paste0("V", seq(nc+1))

	res0 = GenomicDataStream::lmFitResponses(DelayedArray(Y), X)
	res1 = fastLinReg::glmFitResponses(Y, X, "gaussian")
	res2 = GenomicDataStream::glmFitResponses(DelayedArray(Y), X, "gaussian")
	isSame(res0, res1)
	isSame(res1, res2)

	# PROBIT
	for(i in seq(nrow(Y))){
		Y[i,] = sample(2, n, replace=TRUE) - 1
	}
	res1 = fastLinReg::glmFitResponses(Y, X, "probit")
	res2 = GenomicDataStream::glmFitResponses(DelayedArray(Y), X, "probit")
	isSame(res1, res2)

	# LOGIT
	res1 = fastLinReg::glmFitResponses(Y, X, "logit")
	res2 = GenomicDataStream::glmFitResponses(DelayedArray(Y), X, "logit")
	isSame(res1, res2)

	# POISSON
	for(i in seq(nrow(Y))){
		Y[i,] = rnegbin(n, 5, 3) 
	}
	res1 = fastLinReg::glmFitResponses(Y, X, "poisson")
	res2 = GenomicDataStream::glmFitResponses(DelayedArray(Y), X, "poisson")
	isSame(res1, res2)

	# NB
	res1 = fastLinReg::glmFitResponses(Y, X, "nb")
	res2 = GenomicDataStream::glmFitResponses(DelayedArray(Y), X, "nb")
	isSame(res1, res2)


}













