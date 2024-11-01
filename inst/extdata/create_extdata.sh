

# Gabriel Hoffman
# Oct 22, 2024
# Create genotype files in VCF/BCF/BGEN formats

# GenomicDataStream should support .csi index

cd /hpc/users/hoffmg01/.Rlib/R_433/GenomicDataStream/extdata

VCF=/hpc/users/hoffmg01/.Rlib/R_433/BinaryDosage/extdata/set1a.vcf.gz 

ml tabix bcftools parallel #qctool/2.0.6 vcftools/0.1.15
export BGENIX=/hpc/packages/minerva-centos7/bgen/2020-03-13/bin/bgenix

# rm -f test*

# VCF
zcat $VCF | bgzip > test.vcf.gz
tabix -fp vcf test.vcf.gz

# remove phase information
# this gives a BGEN with 3 categories instead of 4
bcftools +setGT test.vcf.gz -- -t a -n u | bgzip > test_noph.vcf.gz
tabix -fp vcf test_noph.vcf.gz

# BCF
bcftools convert test.vcf.gz -o test.bcf
bcftools index test.bcf
bcftools convert test_noph.vcf.gz -o test_noph.bcf
bcftools index test_noph.bcf

# PGEN
#./plink2 --vcf test.vcf.gz dosage=GP-force --make-pgen vzs --out test
plink2 --vcf test.vcf.gz dosage=GP-force --make-pgen --out test
plink2 --vcf test.vcf.gz dosage=GP-force --make-bed --out test


# BGEN with phased data
# reference allele first
./plink2 --vcf test.vcf.gz dosage=GP-force --export bgen-1.1 ref-first --out test_v1.1
./plink2 --vcf test.vcf.gz dosage=GP-force --export bgen-1.2 ref-first bits=8 --out test_v1.2_8bits
./plink2 --vcf test.vcf.gz dosage=GP-force --export bgen-1.3 ref-first bits=8 --out test_v1.3_8bits
./plink2 --vcf test.vcf.gz dosage=GP-force --export bgen-1.2 ref-first bits=16 --out test_v1.2_16bits
./plink2 --vcf test.vcf.gz dosage=GP-force --export bgen-1.3 ref-first bits=16 --out test_v1.3_16bits


# BGEN with UNphased data
# reference allele first
./plink2 --vcf test_noph.vcf.gz dosage=GP-force --export bgen-1.1 ref-first --out test_noph_v1.1
./plink2 --vcf test_noph.vcf.gz dosage=GP-force --export bgen-1.2 ref-first bits=8 --out test_noph_v1.2_8bits
./plink2 --vcf test_noph.vcf.gz dosage=GP-force --export bgen-1.3 ref-first bits=8 --out test_noph_v1.3_8bits
./plink2 --vcf test_noph.vcf.gz dosage=GP-force --export bgen-1.2 ref-first bits=16 --out test_noph_v1.2_16bits
./plink2 --vcf test_noph.vcf.gz dosage=GP-force --export bgen-1.3 ref-first bits=16 --out test_noph_v1.3_16bits

# BGEN index
ls *.bgen | parallel $BGENIX -g {} -index -clobber

# retain 1 .sample file
mv test_v1.1.sample tmp
rm -f *.sample *.log
mv tmp test.sample

# $BGENIX -g test_noph_v1.3_16bits.bgen -list


# scp sklar1:"/hpc/users/hoffmg01/.Rlib/R_433/GenomicDataStream/extdata/test*" .


# test output
# file = "test_v1.2_8bits.bgen"
# res = rbgen::bgen.load( file, ranges = data.frame( chromosome = '1', start = 0, end = 14444), max_entries_per_sample=4 )


# file = "test_noph_v1.2_8bits.bgen"
# res = rbgen::bgen.load( file, ranges = data.frame( chromosome = '1', start = 0, end = 14444), max_entries_per_sample=3 )




# bgenix -g test_v1


# # Write GEN with dosage from DS
# ./plink2 --vcf test.vcf.gz dosage=DS --double-id --export oxford --out test 

# # GEN files can have either 5 or 6 header columns
# # plink writes 5, but qctool expects 6 in order to name chromosomes
# # Otherwise the chromsome entry is blank in the .bgi file
# # and reading with GenomicDataStream fails
# # https://www.cog-genomics.org/plink/2.0/input#oxford
# # Also, I couldn't get qctool to convert from VCF to BGEN directly
# sed -i 's/^1/1 1/g' test.gen 

# qctool -g test.gen -s test.sample -bgen-bits 8 -ofiletype bgen_v1.2 -og test_v1.2_8bits.bgen

# bgenix -g test_v1.2_8bits.bgen -index -clobber
# bgenix -g test_v1.2_8bits.bgen -list


# ./plink2 --vcf test.vcf.gz dosage=GP-force --export bgen-1.2 bits=16 --out test_v1.2_16bits


# bgenix -g test_v1.2_16bits.bgen -index -clobber
# bgenix -g test_v1.2_16bits.bgen -list




# 	file = "test_v1.2_16bits.bgen"
# 	GenomicDataStream:::test_bgen(file, "DS", chunkSize=234, region="1:0-10000")


# 	file = "test_v1.2_16bits.bgen"
# 	D = rbgen::bgen.load( file, ranges = data.frame( chromosome = '1', start = 0, end = 144000044 ), max=4)



# 	D$data[8,1,,drop=FALSE]




# ./plink2 --vcf test.vcf.gz dosage=GP-force --export vcf --out test2

# # 3 categories
# v = c(0, 0.195, 0.805)
# v = c(0.09750515,0.9024949,0.09750515,0.9024949)

# # 0/1:0.897:0.103,0.897,0




# qctool -g test2.vcf.gz -bgen-bits 8 -ofiletype bgen_v1.1 -og test_v1.1_8bits.bgen



# # BCF
# bcftools convert test.vcf.gz -o test.bcf
# bcftools index test.bcf

# # Write GEN with dosage from DS
# ./plink2 --vcf test.vcf.gz dosage=DS --double-id --export oxford --out test 


# qctool -g test.vcf.gz -bgen-bits 8 -ofiletype bgen_v1.1 -og test_v1.1_8bits.bgen



# # write BGEN 8 bits
# qctool -g test.gen -s test.sample -bgen-bits 8 -ofiletype bgen_v1.1 -og test_v1.1_8bits.bgen
# qctool -g test.gen -s test.sample -bgen-bits 8 -ofiletype bgen_v1.2 -og test_v1.2_8bits.bgen

# bgenix -g test_v1.2_8bits.bgen -index
# bgenix -g test_v1.2_8bits.bgen -list


# # write BGEN 16 bits
# qctool -g test.gen -s test.sample -bgen-bits 16 -ofiletype bgen_v1.1 -og test_v1.1_16bits.bgen
# qctool -g test.gen -s test.sample -bgen-bits 16 -ofiletype bgen_v1.2 -og test_v1.2_16bits.bgen

# # index BGEN file
# bgenix -g test_v1.1_8bits.bgen -index -clobber
# bgenix -g test_v1.2_8bits.bgen -index -clobber
# bgenix -g test_v1.1_16bits.bgen -index -clobber
# bgenix -g test_v1.2_16bits.bgen -index -clobber



# bgenix -g test_v1.1_8bits.bgen -vcf | cut -f1




# qctool -g test.gen -s test.sample -bgen-bits 8 -ofiletype bgen_v1.1 -og test_v1.1_8bits.bgen


# zcat test.vcf.gz > test.vcf
# qctool -g test.vcf -ofiletype gen -og test.gen

# qctool -g test.vcf.gz -vcf-genotype-field GP -og converted.bgen


# zcat $VCF | vcf-convert -v 4.2 | bgzip > test.vcf.gz
# tabix -fp vcf test.vcf.gz

 
# qctool -g test.vcf.gz -bgen-bits 8 -ofiletype bgen_v1.1 -og test_v1.1_8bits.bgen





# ./plink2 --vcf test.vcf.gz dosage=GP-force --double-id --export vcf-4.2 --out test2 
# bgzip -f test2.vcf
# tabix -fp vcf test2.vcf.gz
 


# sed -i 's/^1/chr9/g' test.vcf
# sed -i 's/<ID=9>/<ID=chr9>/g' test.vcf

# # VCF
# gunzip -f tmp_test.vcf.gz

# cat tmp_test.vcf | sed 's/^1/21/g' > tmp_test2.vcf

# bgzip tmp_test.vcf
# tabix -p vcf tmp_test.vcf.gz
# bgzip tmp_test2.vcf
# tabix -p vcf tmp_test2.vcf.gz


# bcftools concat tmp_test.vcf.gz tmp_test2.vcf.gz | bgzip > test.vcf.gz
# tabix -p vcf test.vcf.gz


# OLD
# # write BGEN with 8 bits with dosage from DS
# plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.1 --out test_v1.1
# plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.2 bits=8 --out test_v1.2_8bits
# plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.3 bits=8 --out test_v1.3_8bits

# # write BGEN with 16 bits with dosage from DS
# plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.2 bits=16 --out test_v1.2_16bits
# plink2 --vcf test.vcf.gz dosage=DS --export bgen-1.3 bits=16 --out test_v1.3_16bits

# rm -f *.log test_v*.sample

# alias bgenix=/hpc/packages/minerva-centos7/bgen/2020-03-13/bin/bgenix
# bgenix -g test_v1.1.bgen -index
# bgenix -g test_v1.2_8bits.bgen -index
# bgenix -g test_v1.3_8bits.bgen -index
# bgenix -g test_v1.2_16bits.bgen -index
# bgenix -g test_v1.3_16bits.bgen -index
