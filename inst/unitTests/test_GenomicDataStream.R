


test_sample_order = function(){

	library(GenomicDataStream)
	library(RUnit)

	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
	files = list.files(dirname(file), "(vcf.gz|bcf|bgen|pgen)$", full.names=TRUE)

	# remove BGEN v1.1, since this doesn't store sample names
	files = grep("1.1.bgen", files, value=TRUE, invert=TRUE)

	# read full data
	gds <- GenomicDataStream(file, "GT", initialize = TRUE)
	dat <- getNextChunk(gds)
	ids = rownames(dat$X)

	# check sample ordering based on sample query
	# BGEN: uses order from query
	# VCF/BCF/PGEN: uses order in file
	res = c()
	for(file in files){
	# 	cat(basename(file), ": ")	
		ids.rnd = sample(ids, length(ids), replace=FALSE)

		gds1 <- GenomicDataStream(file, "GT", initialize = TRUE, samples=ids.rnd)
		dat1 <- getNextChunk(gds1)
		ids1 = rownames(dat1$X)

		allSame = sum(ids == ids1) == 60
		res = c(res, allSame)
	}
	names(res) = basename(files)

	# BGEN not equal
	a = all(!res[grep("bgen$", names(res))])
	checkEquals(a, TRUE)

	# VCF/BCF/PGEN are equal
	a = all(res[grep("(vcf.gz|bcf|pgen)$", names(res))])
	checkEquals(a, TRUE)


	for(file in files){
		# cat(basename(file), ": ")	
		ids.rnd = sample(ids, length(ids), replace=FALSE)

		gds1 <- GenomicDataStream(file, "GT", initialize = TRUE, samples=ids.rnd)
		dat1 <- getNextChunk(gds1)

		# reorder using original order
		ids1 = rownames(dat1$X[ids,])

		checkEquals(ids, ids1)
	}

	# check that getSampleNames() and getNextChunk()
	# give same sample order
	ids.rnd = sample(ids, length(ids), replace=FALSE)

	for(file in files){
		# cat(basename(file), ": ")	

		gds1 <- GenomicDataStream(file, "GT", initialize = TRUE, samples=ids.rnd)
		dat1 <- getNextChunk(gds1)

		checkEquals(getSampleNames(gds1), rownames(dat1$X))
	}

}





test_multiple_GenomicDataStream = function(){

	library(GenomicDataStream)
	library(RUnit)

	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
	files = list.files(dirname(file), "(vcf.gz|bcf|bgen|pgen)$", full.names=TRUE)
	files = grep("1.1.bgen", files, value=TRUE, invert=TRUE)

	for(file in files){
		# cat(file, "\n")
		# read full file
		obj <- GenomicDataStream(file, field = "GT", chunkSize = 30, init=TRUE)
		dat1 <- getNextChunk(obj)  

		# read subset without re-initializing
		gds = setRegion(obj, "1:15000-17000")
		dat2 <- getNextChunk(gds)   

		# equal to first set
		checkEquals(c(dat2$info), c(dat1$info[6:8,]))
		checkEquals(dat2$X, dat1$X[,6:8])

		# enpty now
		checkEquals(length(getNextChunk(gds)), 0)

		# read subset without re-initializing
		gds = setRegion(obj, "1:11000-14001")
		dat2 <- getNextChunk(gds) 

		# equal to first set
		checkEquals(c(dat2$info), c(dat1$info[2:5,]))
		checkEquals(dat2$X, dat1$X[,2:5])

		# enpty now
		checkEquals(length(getNextChunk(gds)), 0)
	}
}


test_minVariance_filter = function(){

	library(GenomicDataStream)
	library(RUnit)

	file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")
	files = list.files(dirname(file), "(vcf.gz|bcf|bgen|pgen)$", full.names=TRUE)
	MAF = .3
	minVar = 2*MAF*(1-MAF)

	for(file in files){
		# cat(file, "\n")				

		gds1 <- GenomicDataStream(file, "GT", init=TRUE)
		dat1 <- getNextChunk(gds1)
		i = (apply(dat1$X, 2, var) > minVar)
		which(i)

		gds2 <- GenomicDataStream(file, "GT", init=TRUE, MAF=MAF)
		dat2 <- getNextChunk(gds2)

		checkEquals(dat1$X[,i], dat2$X)
		checkEquals(c(dat1$info[i,]), c(dat2$info))
	}
}




test_DataTable = function(){
	# q()
	# R
	suppressPackageStartupMessages({
	library(GenomicDataStream)
	})

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")

	# file <- system.file("extdata", "test.pvar", package = "GenomicDataStream")
	# res1 = GenomicDataStream:::test_DataTable( file, "#CHROM")


	# file <- system.file("extdata", "test.psam", package = "GenomicDataStream")
	# res2 = GenomicDataStream:::test_DataTable( file, "#IID")


	# file <- system.file("extdata", "test.fam", package = "GenomicDataStream")
	# res3 = GenomicDataStream:::test_DataTable( file, "")


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

	# minVariance = NaN means no filtering
	gdsObj = GenomicDataStream(file, region=reg, samples=ids_str,chunkSize=3, initialize=TRUE, minVariance=NaN)
	dat = getNextChunk(gdsObj)

	pg = NewPgen(file, raw_sample_ct=60, sample_subset=idx)
	X = ReadList(pg, 2:3)

	value = max(abs(X - dat$X), na.rm=TRUE)

	checkEqualsNumeric(value, 0)
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
		
		# print(dat$info)
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
		# cat(file, "\n")
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
	# fit2 = GenomicDataStream:::test_lm(X,y)
	# checkEqualsNumeric(coef(fit), fit2$coefficients)
	
	# Analysis in GenomicDataStream
	#-------------------------------

	files = list.files(dirname(file), "(vcf.gz|bcf|bgen)$", full.names=TRUE) # |pgen

	X_design = matrix(1, ncol(X_all))
	w = y
	w[] = 1

	# All regions
	###############
	resList = lapply(files, function(file){

		# cat(file, "\n")
		# rm(dat, res)

		# test dosages
		gds = GenomicDataStream(file, "DS", chunkSize=10000, initialize=TRUE, missingToMean=TRUE)
	    dat = getNextChunk(gds)

		# dat$X[1:3, 1:3]
		# t(X_all[1:3,1:3])

		checkEqualsNumeric(t(dat$X), X_all, tol=1e-4) 
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

		# cat(file, "\n")

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
		})

}

