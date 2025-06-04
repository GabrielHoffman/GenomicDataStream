
q()
R
library(GenomicDataStream)

file = "~/Downloads/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

region = "chr21:0-18030588"
gds = GenomicDataStream(file, "GT", chunkSize = 1000, region=region, init=TRUE)

system.time(res <- PCAstream(gds, k=2))

835.561

# barplot(res$timing)


gds = GenomicDataStream(file, "GT", chunkSize = 1e7, region=region, init=TRUE)

system.time(dat <- getNextChunk(gds))
standardize_in_place(dat$X)

dim(dat$X)

k = 2

# system.time(s0 <- svd(dat$X, k, k))
system.time(s1 <- pcaone::pcaone(dat$X, k = k))
system.time(s2 <- RSpectra::svds(dat$X, k = k))
system.time(s3 <- irlba::irlba(dat$X, k = k))

res$d
s1$d

par(mfrow=c(1,2))
plot(res$v)
plot(s1$u)



q()
R
library(GenomicDataStream)

file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

gds <- GenomicDataStream(file, "DS", chunkSize = 3, init=TRUE)

res <- PCAstream(gds, k=10)

res



