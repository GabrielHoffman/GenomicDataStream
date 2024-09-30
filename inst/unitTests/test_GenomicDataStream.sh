
test_vcfstream = function(){

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

	res = lucida:::extractVcf( file, "DS", "." )
	checkEqualsNumeric(t(res$X), X_all, tol=1e-7)

	res2 = vcfppR::vcftable(file, ".", format="DS")
	checkEqualsNumeric(res2$DS, X_all, tol=1e-7)

	res3 = lucida:::extractVcf_chunks( file, "DS", "." )
	checkEqualsNumeric(res$X[,3:4], res3$X, tol=1e-7)

	# range
	#-------------

	# VariantAnnotation
	region = "1:11000-13000"
	gr = GRanges(1, IRanges(11000, 13000))
	vcf <- suppressWarnings(readVcf(file, "hg19", param=gr))
	X_all = geno(vcf)[["DS"]]

	res = lucida:::extractVcf( file, "DS", region )
	checkEqualsNumeric(t(res$X), X_all, tol=1e-7)

	res3 = lucida:::extractVcf_chunks( file, "DS", region )
	checkEqualsNumeric(res$X[,3], res3$X, tol=1e-7)

	# subset samples
	#---------------
	idx = c(6, 1, 9, 4)
	sampleIds = paste0(rownames(res$X)[idx], collapse=",")
	res2 = lucida:::extractVcf( file, "DS", region, sampleIds )
	checkEqualsNumeric(res2$X, res$X[sort(idx),], tol=1e-7)

	res3 = lucida:::extractVcf_chunks( file, "DS", region, sampleIds  )
	checkEqualsNumeric(res$X[sort(idx),3], res3$X, tol=1e-7)


	# DP is an integer
	##################
	file <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
    
    field = "DP"
	vcf <- suppressWarnings(readVcf(file))
	X_all = geno(vcf)[[field]]
	res = lucida:::extractVcf( file, field, "." )
	checkEqualsNumeric(t(res$X), X_all)

	res2 = lucida:::extractVcf_chunks( file, field, "." )
	checkEqualsNumeric(t(res2$X), X_all[3:4,])

	# GT must be converted to integer
	#################################
	field = "GT"	

	region = "chr21:5030082-5030105"
	res = lucida:::extractVcf( file, field, region)
	res2 = lucida:::extractVcf_chunks( file, field, region)

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
	res2 = lucida:::extractVcf( file, field, region, missingToMean=TRUE)
	a = colMeans(res$X, na.rm=T)
	b = colMeans(res2$X)
	checkEqualsNumeric(a,b)

	res3 = lucida:::extractVcf_chunks( file, field, region, missingToMean=TRUE)
	a = colMeans(res$X, na.rm=T)
	b = colMeans(res3$X)
	checkEqualsNumeric(a[3], b)



	# Check errors
	#-------------

	# res = lucida:::extractVcf( file, "GT", "." )
	# Error: GT is not supported for multi-allelic site
	# chr21:5030319 chr21:5030319:C:G,T C G,T

	# res = lucida:::extractVcf( file, "GT", "343" )
	# Error: region was not found! make sure the region format is correct

     
}




