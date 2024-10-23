

# Gabriel Hoffman
# Oct 22, 2024
# Create genotype files in VCF/BCF/BGEN formats

# GenomicDataStream should support .csi index

cd /Users/gabrielhoffman/workspace/repos/GenomicDataStream/inst/extdata

\cp -f /Users/gabrielhoffman/prog/R-4.4.0/library/BinaryDosage/extdata/set1a.vcf.gz test.vcf.gz

# VCF
gunzip test.vcf.gz
bgzip test.vcf
tabix -p vcf test.vcf.gz

# BCF
bcftools convert test.vcf.gz -o test.bcf

# index BCF
bcftools index test.bcf

# Write GEN with dosage from DS
plink2 --vcf test.vcf.gz dosage=DS --export oxford --out test

# write BGEN with 8 bits with dosage from DS
plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.1 --out test_v1.1
plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.2 bits=8 --out test_v1.2_8bits
plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.3 bits=8 --out test_v1.3_8bits

# write BGEN with 16 bits with dosage from DS
plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.2 bits=16 --out test_v1.2_16bits
plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.3 bits=16 --out test_v1.3_16bits

rm -f *.log test_v*.sample

# alias bgenix=/hpc/packages/minerva-centos7/bgen/2020-03-13/bin/bgenix
# bgenix -g test_v1.1 -index
# bgenix -g test_v1.2_8bits -index
# bgenix -g test_v1.3_8bits -index
# bgenix -g test_v1.2_16bits -index
# bgenix -g test_v1.3_16bits -index







