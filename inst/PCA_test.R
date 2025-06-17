


ml tabix bcftools plink2/v2.00a3.3
cd /sc/arion/scratch/hoffmg01

remotes::install_github("Zilong-Li/PCAoneR")

# Get data
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

# Filter by MAF just to get a smaller file to play with
# FILTER="(AF > 0.10) & (AF < 0.90) & (TYPE='snp') & (ALT='A')"
FILTER="(AF > 0.01) & (AF < 0.99) & (TYPE='snp') & (ALT='A')"

bcftools view -i "(AF > 0.01) & (AF < 0.99) & (TYPE='snp')" 1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | bcftools view -O b -o 1kg_chr21.bcf
bcftools index 1kg_chr21.bcf

plink2 --bcf 1kg_chr21.bcf --allow-extra-chr --make-bed --out 1kg_chr21

plink2 --bcf 1kg_chr21.bcf --allow-extra-chr --make-pgen --out 1kg_chr21





# Combine all of 1KG
ml dataark parallel
PARENT=/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/VCF
FILES=$(seq 1 22 | parallel ls $PARENT/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz)

# keep common biallelic SNPs
bcftools concat $FILES | bcftools view -i "(AF > 0.01) & (AF < 0.99) & (TYPE='snp')" -m 2 -M 2 | bcftools view -O b -o 1kg.bcf


plink --bcf 1kg.bcf --allow-extra-chr --maf 0.1 --make-bed --out 1kg


plink2 --bcf 1kg_chr21.bcf --allow-extra-chr --make-bed --out 1kg_chr21
plink2 --bcf 1kg_chr21.bcf --allow-extra-chr --make-pgen --out 1kg_chr21


plink2 --bcf 1kg_chr1.bcf --allow-extra-chr --make-bed --out 1kg_chr1
plink2 --bcf 1kg_chr1.bcf --allow-extra-chr --make-pgen --out 1kg_chr1

# Flash PCA scales poorly with number of PCs
plink2 --bcf 1kg_chr1.bcf --pca 20 approx



# VCF to filtered BCF and BED
FILE=/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
bcftools view -i "(AF > 0.01) & (AF < 0.99) & (TYPE='snp')" -m 2 -M 2 $FILE | bcftools view -O b -o 1kg_chr1.bcf
bcftools index 1kg_chr1.bcf
plink --bcf 1kg_chr1.bcf --allow-extra-chr --make-bed --out 1kg_chr1



cd ~/Downloads



library(GenomicDataStream)
file = "test.vcf.gz"
gds = GenomicDataStream(file, field="GT", init=TRUE)

library(vcfppR)
res = vcftable(file)



q()
R

suppressPackageStartupMessages({
library(GenomicDataStream)
library(DelayedArray)
library(beachmat)
library(RUnit)
})

n <- 10000
m <- 25
k = 3
set.seed(1)
Y <- matrix(rpois(n * m, 100), m, n)
Y <- DelayedArray(Y)

dcmp = svd(scale(Y) / sqrt(nrow(Y)-1), k, k)

res = PCAstream(Y, k=k, 5000, verbose=TRUE)


checkEqualsNumeric(dcmp$d[seq(k)], res$d)
checkEqualsNumeric(abs(dcmp$u), abs(res$u))
checkEqualsNumeric(abs(dcmp$v), abs(res$v))



perfMetric(dcmp$u, res$u, k, "MEV")
perfMetric(dcmp$u, res$u, k, "minSSE")









q()
R

suppressPackageStartupMessages({
library(GenomicDataStream)
library(DelayedArray)
library(muscat)
library(SingleCellExperiment)
library(RUnit)
})

data(example_sce)

k = 50

Y = logcounts(example_sce) + 23

dcmp = svd(scale(Y), k, k)

res = PCAstream(t(Y), k=k, 1110, verbose=TRUE, p=12)

checkEqualsNumeric(dcmp$d[seq(k)], res$d)
checkEqualsNumeric(abs(dcmp$u), abs(res$u))
checkEqualsNumeric(abs(dcmp$v), abs(res$v))


perfMetric(dcmp$u, res$u, k, "MEV")
perfMetric(dcmp$u, res$u, k, "minSSE")


values = sapply(seq(k), function(i){
	perfMetric(dcmp$u, res$u, i, "MEV")
	})
plot(seq(k), values)




file = "/Users/gabrielhoffman/prog/R-4.4.2/library/GenomicDataStream/extdata/test_noph_v1.3_16bits.bgen"
gds <- GenomicDataStream::GenomicDataStream(file, init=TRUE)



q()
R

suppressPackageStartupMessages({
library(GenomicDataStream)
library(DelayedArray)
library(muscat)
library(SingleCellExperiment)
library(RUnit)
library(ExperimentHub)
library(dreamlet)
library(scater)
})

# Download data, specifying EH2259 for the Kang, et al study
eh <- ExperimentHub()
sce <- eh[["EH2259"]]

sce$ind = as.character(sce$ind)

# only keep singlet cells with sufficient reads
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
sce <- sce[, colData(sce)$multiplets == "singlet"]

# compute QC metrics
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]

# set variable indicating stimulated (stim) or control (ctrl)
sce$StimStatus <- sce$stim

prior.count = 1
lib.size = colSums2(counts(sce))
logcounts(sce) = t(log2(t(counts(sce) + prior.count)) - log2(lib.size) + log2(1e6))


Y = logcounts(sce)

k = 20

# dcmp <- svd(scale(Y), k, k)
dcmp <- irlba::irlba(scale(Y), k, k)

res = PCAstream(t(Y), k=k, 1000, verbose=TRUE, p=10)

checkEqualsNumeric(dcmp$d[seq(k)], res$d)
checkEqualsNumeric(abs(dcmp$u), abs(res$u))
checkEqualsNumeric(abs(dcmp$v), abs(res$v))

times



q()
R
suppressPackageStartupMessages({
library(DelayedArray)
library(zellkonverter)
library(SingleCellExperiment)
library(RUnit)
library(DelayedMatrixStats)
library(BiocSingular)
library(GenomicDataStream)
library(dreamlet)
})

file = "/sc/arion/projects/CommonMind/leed62/GENESIS/GEN_A2/GEN_A2_pass3_anno.h5ad"
sce_in = readH5AD(file, use_hdf5=TRUE, verbose=TRUE, raw=TRUE, uns=FALSE, obsp=FALSE, obsm=FALSE)
sce <- swapAltExp(sce_in, "raw") # use raw as main
counts(sce) <- assay(sce, "X") 

 
sce$total_counts = colSums2(counts(sce))
lib.size = sce$total_counts
prior.count = 1
logcounts(sce) <- t(log2(t(counts(sce) + prior.count)) - log2(lib.size) + log2(1e6))

k = 40
system({
res1 <- PCAstream(t(logcounts(sce)), k=k, chunkSize=1000, threads=6)
})

q()
R
suppressPackageStartupMessages({
library(DelayedArray)
library(zellkonverter)
library(SingleCellExperiment)
library(RUnit)
library(DelayedMatrixStats)
library(BiocSingular)
library(GenomicDataStream)
})

# devtools::reload("/hpc/users/hoffmg01/.Rlib/R_441_bioc320/GenomicDataStream")

file = "~/Downloads/4e6932db-5a78-40e4-b961-f87f66ba139a.h5ad"
# file = "~/4e6932db-5a78-40e4-b961-f87f66ba139a.h5ad"
sce_in = readH5AD(file, use_hdf5=TRUE, raw=TRUE, verbose=TRUE, uns=FALSE, obsp=FALSE, obsm=FALSE)
sce <- swapAltExp(sce_in, "raw") 
counts(sce) <- assay(sce, "X") 
sce$total_counts = colSums2(counts(sce))

lib.size = sce$total_counts
prior.count = 1
logcounts(sce) <- t(log2(t(counts(sce) + prior.count)) - log2(lib.size) + log2(1e6))

sce = sce[, 1:6000]

k = 4
res1 <- PCAstream(t(logcounts(sce)), k=k, chunkSize=1000, threads=6)
res2 <- PCAstream(sce, assay="logcounts", k=k, chunkSize=1000, threads=6)



X = t(logcounts(sce))
X_scale = scale(X) / sqrt(nrow(X) -1)
X_scale = as.matrix(X_scale)

system.time(
	dcmp <- irlba::irlba(X_scale, k, k)
	)


checkEqualsNumeric(dcmp$d[seq(k)], res1$d, tol=1e-3)
checkEqualsNumeric(abs(dcmp$u), abs(res1$u), tol=1e-3)
checkEqualsNumeric(abs(dcmp$v), abs(res1$v), tol=1e-3)






res1 <- PCAstream(sce2, k=10, chunkSize=500, threads=1)

i = 1:5000
a = as.matrix(logcounts(sce)[,i])
a = as.matrix(logcounts(sce)[i,])



rs = rowSums2(counts(sce2))
keep = (rs > 200)


k = 50
X_scale = scale(as.matrix(logcounts(sce2)))
system.time(
	dcmp <- irlba::irlba(X_scale, k, k)
	)


res1 <- PCAstream(sce[keep, seq(10000)], k=10, chunkSize=5000, threads=1)



system.time({
	resone <- pcaone::pcaone(X_scale, k=k)
})

system.time({
	cs.out <- runIrlbaSVD(Y, k=k, scale=colSds(Y), 
    center=colMeans2(Y))
})


system.time(	
	res1 <- PCAstream(sce, k=k, assay="sdf")
)

system.time(	
	res2 <- PCAstream(Y, k=k)
)


times

checkEqualsNumeric(dcmp$d[seq(k)], res$d)
checkEqualsNumeric(abs(dcmp$u), abs(res$u))
checkEqualsNumeric(abs(dcmp$v), abs(res$v))


par(mfrow=c(1,2))
plot(dcmp$v[,1:2])
plot(res$v[,1:2])


q()
R
library(muscat)
library(SingleCellExperiment)
library(GenomicDataStream)

data(example_sce)
sce <- example_sce

# Normalize expression with log2 counts per million
prior.count <- 1
lib.size <- colSums2(counts(sce))
logcounts(sce) <- t(log2(t(counts(sce) + prior.count)) - log2(lib.size) + log2(1e6))

# PCA on SingleCellExperiment with assay logcounts
res <- PCAstream(sce, k=5, assay="logcounts")


  trace("PCAstream", browser, exit=browser, signature = c("SingleCellExperiment")) 


FILE=/sc/arion/scratch/hoffmg01/1kg
time ~/build2/PCAone/PCAone --bfile $FILE -m 20



FILE=/sc/arion/scratch/hoffmg01/1kg_chr1
~/build2/PCAone/PCAone --bfile $FILE -m 20


q()
R
suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "~/Downloads/1kg_chr21.bcf"
file = "1kg_chr21.bed"
k = 50
MAF = .48

gds = GenomicDataStream(file, field="GT", init=TRUE, chunkSize=100000, MAF=MAF)

X = c()
while(1){
	dat = getNextChunk(gds)	
	X = cbind(X, dat$X)
  if (atEndOfStream(gds)) break
}

X_scale = scale(X) / sqrt(nrow(X)-1)
dcmp = irlba::irlba(X_scale, k, k)

# res = pcaone::pcaone( X_scale, k=k, shuffle=TRUE )

# Streaming winSVD
# small chunkSize
gds = GenomicDataStream(file, field="GT", init=TRUE, chunkSize=10000, MAF=MAF)

system.time(
res <- PCAstream(gds, k)
)
perfMetric(dcmp$u, res$u, k, "MEV")


values = sapply(seq(k), function(i){
	perfMetric(dcmp$u, res$u, i, "MEV")
	})
plot(seq(k), values, ylim=c(0, 1))


perfMetric(dcmp$u, res$u, k, "minSSE")





trace("pcaone", browser, exit=browser, signature = c("matrix")) 




trace("PCAstream", browser, exit=browser, signature = c("GenomicDataStream")) 



trace("PCAstream", browser, exit=browser, signature = c("ANY")) 



q()
R




suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "~/Downloads/1kg_chr1.bed"
k = 100
gds = GenomicDataStream(file, "GT", chunkSize=100000, MAF=0.05 )

system.time({
res <- PCAstream(gds, k)
})





res <- PCAstream(gds, k, verbose=FALSE, shuffle=TRUE)




sObj = summaryChunks(gds)

gds2 = GenomicDataStream(file, "GT", chunkSize=1000, region=paste(sObj$regions, collapse = "," ), init=TRUE)
system.time(
res <- PCAstream(gds, k, verbose=TRUE, shuffle=TRUE )
)



system.time(
res <- PCAstream(gds, k, verbose=TRUE)
)



plot(res2$u[,1:2])



q()
R
suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "~/Downloads/1kg_chr1.bed"
gds = GenomicDataStream(file, chunkSize=100000, init=TRUE)






system.time(
a <- summaryChunks(gds)
)


q()
R
suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "~/Downloads/1kg_chr21.pgen"
gds = GenomicDataStream(file, "GT", init=FALSE, chunkSize=10000)

system.time({
res = PCAstream(gds, 40, verbose=TRUE)
})





q()
R
suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "~/Downloads/1kg_chr21.pgen"
gds = GenomicDataStream(file, "GT", init=FALSE, chunkSize=100000)
system.time({
a <- summaryChunks(gds)
})



q()
R
system.time({
GenomicDataStream:::test_DataTable2()
})





  trace("PCAstream", browser, exit=browser, signature = c("GenomicDataStream")) 


suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "~/Downloads/1kg_chr21.bed"
gds = GenomicDataStream(file, init=FALSE, chunkSize=10020)
a <- summaryChunks(gds)


suppressPackageStartupMessages({
library(GenomicDataStream)
})
ids = read.table("~/Downloads/1kg_chr21.fam")$V2
file = "~/Downloads/1kg_chr21.pgen"
gds = GenomicDataStream(file, init=FALSE, chunkSize=1332, samples=paste0(ids[1:2000], collapse=','))
a <- summaryChunks(gds)

repeat {
	a <- summaryChunks(gds)
}


Rcpp::sourceCpp(code = "
#include <RcppArmadillo.h>
using namespace Rcpp ;
// [[Rcpp::depends(RcppArmadillo)]]

class Obj {
	public:
	Obj() {}
	Obj(arma::mat x): x(x) {}

	arma::mat x;
};

void f(Obj &x){
	std::vector<double> matDosage(100);
	arma::mat M(matDosage.data(), 10, 10, false, true);
	x = Obj(M);
}

// [[Rcpp::export]]
List getData(){
	Obj d;
	f(d);
	arma::mat M = d.x;
	List lst = List::create(Named(\"X\") = wrap(M));
	return(lst);
}
")

repeat { 
	cat("a")
	a = rnorm(10333000)
	res = getData()
}


countChunks(gds)



library(pgenlibr)
pvar_path = "1kg_chr21.pvar"
pgen_path = "1kg_chr21.pgen"
pvar <- pgenlibr::NewPvar(pvar_path)
pgen <- pgenlibr::NewPgen(pgen_path, pvar=pvar)

variant_subset = sapply(seq(GetVariantCt(pgen)), 
	function(i){
	GetVariantId(pvar, i)
	})

res = ReadList(pgen, seq(GetVariantCt(pgen)))





q()
R
suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "~/Downloads/1kg_chr21.bed"
gds = GenomicDataStream(file, init=TRUE)
sdf = rnorm(10000)
res = GenomicDataStream:::summarizeChunks_rcpp(gds@ptr)

a <- summaryChunks(gds)


q()
R
suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "~/Downloads/1kg_chr21.bcf"
gds = GenomicDataStream(file, "GT", init=TRUE)
res = GenomicDataStream:::summarizeChunks_rcpp(gds@ptr)


q()
R
suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "/Users/gabrielhoffman/Downloads/1kg_chr21.bcf"
GenomicDataStream:::create_test( file, "GT")


file = "/Users/gabrielhoffman/Downloads/1kg_chr21.bed"
GenomicDataStream:::create_test( file )


gds = GenomicDataStream(file, init=TRUE)
res = GenomicDataStream:::summarizeChunks_rcpp(gds@ptr)


file = "/Users/gabrielhoffman/Downloads/1kg_chr21.bed"
GenomicDataStream:::create_test( file )


GenomicDataStream:::create_test2()


count <- 0
while(1){
	dat = getNextChunk(gds)	
	count <- count + ncol(dat$X)
	cat(paste0("count: ", count, "\n"))
    if (atEndOfStream(gds)) break
}



suppressPackageStartupMessages({
library(GenomicDataStream)
})
file = "~/Downloads/1kg_chr21.bed"
gds = GenomicDataStream(file, init=FALSE, chunkSize=100020)
a <- summaryChunks(gds)










for(i in 1:100){
gds = GenomicDataStream(file, init=TRUE, region="21:0-524126500", chunkSize=i)
a <- summaryChunks(gds)
}

gds = GenomicDataStream(file, init=TRUE)
a <- summaryChunks(gds)


q()
R

suppressPackageStartupMessages({
library(GenomicDataStream)
})
file <- system.file("extdata", "test.bed", package = "GenomicDataStream")
for(i in 1:100){
gds = GenomicDataStream(file, init=TRUE, chunkSize=i)
a <- summaryChunks(gds)
}










system.time(
res2 <- PCAstream(gds, k, verbose=TRUE)
)


 trace("PCAstream", browser, exit=browser, signature = c("GenomicDataStream")) 


dcmp1 = dcmp
dcmp2 = res2

concordance <- function(dcmp1, dcmp2){

	n = nrow(dcmp1$u)
	p = nrow(dcmp1$v)

	# frac = fraction of variance in the ith PC
	# f = sqrt(frac)
	dcmp1$f = with(dcmp1, sqrt(d^2 / p))
	dcmp2$f = with(dcmp2, sqrt(d^2 / p))

	A = with(dcmp1, sqrt(n-1) * u %*% diag(f[1:ncol(u)]))
	B = with(dcmp2, sqrt(n-1) * u %*% diag(f[1:ncol(u)]))

	a = apply(A,2,var)
	apply(B,2,var)

	diag(cov(A,B))
	diag(cor(A,B))

	a = diag(cor(A,B))^2
	b = diag(cov(A,B))^2

	weighted.mean(a,b)


	# res_cca = cancor(A, B)

	# sqrt(mean(res_cca$cor^2))

	# a = A %*% normalize(res_cca$xcoef)
	# b = B %*% normalize(res_cca$ycoef)
	# diag(cor(a,b))
}


normalize = function(H){

	# a = x / sqrt(sum(x^2))
	# sum(a^2)
	apply(H, 2, function(x) x / sqrt(sum(x^2)))
}




# Plot
par(mfrow=c(2,2))
plot(dcmp$d)
plot(dcmp$u[,1:2])
ylim = range(c(res$d, res2$d))
plot(dcmp$d[seq(k)], res$d, xlab="Standard SVD", ylab="Window SVD", col="blue", pch=15, ylim=ylim)
points(dcmp$d[seq(k)], res2$d, col="green3", lwd=3)
legend("topleft", legend = c("winSVD", "winSVDstream"), fill=c("blue", "green3"))
abline(0, 1, col="red")

plot(diag(abs(cor(dcmp$u, res$u))), ylab="Correlation between PC estimates", col="blue", pch=15, ylim=c(0,1))
points(diag(abs(cor(dcmp$u, res2$u))), col="green3", lwd=3 )  
legend("bottomright", legend = c("winSVD", "winSVDstream"), fill=c("blue", "green3"))

concordance(dcmp, res)
concordance(dcmp, res2)

A = dcmp$u %*% diag(dcmp$d[1:ncol(dcmp$u)])
B = res$u %*% diag(res$d)
sqrt(mean(cancor(A, B)$cor^2))


C = res2$u %*% diag(res2$d)
sqrt(mean(cancor(A, C)$cor^2))






cramerw = function(corr, v){

	w <- v / sum(v)
	# sqrt(mean(cancor(dcmp$u, res$v)$cor^2))
	sqrt(weighted.mean(corr^2, w))
}

cramerw(cancor(dcmp$u, res$u)$cor, dcmp$d^2)

A = with(dcmp, u %*% diag(d[1:ncol(u)]))
B = with(res, v %*% diag(d[1:ncol(v)]))
sqrt(mean(cancor(A, B)$cor^2))

B = with(res2, v %*% diag(d[1:ncol(v)]))
sqrt(mean(cancor(A, B)$cor^2))

concordance(dcmp, res)
concordance(dcmp, res2)




devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")






# bcftools view -O b -o panel_full.bcf panel_full.vcf.gz

q()
R
library(GenomicDataStream)
file = "~/Downloads/panel_full.bcf"

region = "22"
MAF = 0.0
minVariance = 2 * (1 - MAF) * MAF
k = 10
gds = GenomicDataStream(file, field="GT", init=TRUE, region=region, chunkSize=100000, minVariance = minVariance)
dat = getNextChunk(gds)
res = winSVD( scale(dat$X), k=k, B=4 )
dcmp = svd(scale(dat$X), k, k)

gds = GenomicDataStream(file, field="GT", init=TRUE, chunkSize=100, region=region, minVariance = minVariance)

res2 = winSVDstream(gds, k, B=4, p=1)


par(mfrow=c(1,2))
plot(dcmp$d[seq(k)], res$d)
abline(0, 1, col="red")

plot(dcmp$d[seq(k)], res2$d)
abline(0, 1, col="red")



# random order
locs = getVariantLocations(gds)
locs2 = sample(locs, length(locs), replace=FALSE)
region = paste(locs2, collapse=",")

# evenly spaced
# locs = getVariantLocations(gds)
# locs2 = locs[seq(1,length(locs), length.out=15000)]
# region = paste(locs2, collapse=",")


# reading with region is slow
# see reader->getStatus( region )
# But if chunks are thinned, does it need to be permuated
gds2 = GenomicDataStream(file, field="GT", init=TRUE, region=region, chunkSize=5000)

locs3 = getVariantLocations(gds2)



res2 = winSVDstream(gds2, 5, B=4)

# TODO
# randomize features
# handle fixed 



devtools::reload("/Users/gabrielhoffman/workspace/repos/GenomicDataStream")


# slow region access
q()
R
library(GenomicDataStream)
file = "~/Downloads/panel_full.bcf"

region = "22"

gds = GenomicDataStream(file, field="GT", init=TRUE, chunkSize=10000, region=region)

locs = getVariantLocations(gds)
region = paste(locs, collapse=",") # many chunks
# region = "22:10534570-50807822" # one chunk 

gds2 = GenomicDataStream(file, field="GT", init=TRUE, region=region, chunkSize=5000)

locs3 = getVariantLocations(gds2)

obj = getNextChunk(gds2)





