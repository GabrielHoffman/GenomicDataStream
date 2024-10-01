/***********************************************************************
 * @file        export.cpp
 * @author      Gabriel Hoffman
 * @email       gabriel.hoffman@mssm.edu
 * @brief       Expose GenomicDataStream library to R
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifdef ARMA
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifdef EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 

#include <vcfstream.h>
#include <bgenstream.h>
#include <vcfpp.h>
#include <DelayedStream.h>

using namespace std;
using namespace vcfpp;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::vec test_GT( const std::string &vcffile,
			const std::string &region){
	// "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
	BcfReader vcf( vcffile, region );
    BcfRecord var(vcf.header); // construct a variant record

    vector<char> gt; // genotype can be bool, char or int type
    vector<int> hetsum(vcf.nsamples, 0);

    vec v;

    while (vcf.getNextVariant(var)) {

        var.getGenotypes(gt);

        v = conv_to< vec >::from(gt);


        if (!var.isSNP() || !var.isNoneMissing()) continue; 

        assert(var.ploidy()==2); // make sure it is diploidy


    }

    return v;

}


// [[Rcpp::export]]
arma::vec test_DS( const std::string &vcffile,
			const std::string &region){
	// "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
	BcfReader vcf( vcffile, region );
    BcfRecord var(vcf.header); // construct a variant record

    vector<float> gt; // genotype can be bool, char or int type
    vector<int> hetsum(vcf.nsamples, 0);

    vec v;

    while (vcf.getNextVariant(var)) {

        var.getFORMAT("GT", gt);

        v = conv_to< vec >::from(gt);


        if (!var.isSNP() || !var.isNoneMissing()) continue; 

        assert(var.ploidy()==2); // make sure it is diploidy

    }

    return v;

}

// convert to List
DataFrame toDF( const VariantInfo &vInfo){

    // return created data frame
    return DataFrame::create(
                    Named("CHROM") = Rcpp::wrap(vInfo.CHROM),
                    Named("POS") = Rcpp::wrap(vInfo.POS),
                    Named("ID") = Rcpp::wrap(vInfo.ID),
                    Named("REF") = Rcpp::wrap(vInfo.REF),
                    Named("ALT") = Rcpp::wrap(vInfo.ALT),
                    _["stringsAsFactors"] = false);
}


// [[Rcpp::export]]
List extractVcf( 
            const std::string &file,
            const std::string &field,
            const std::string &region = "",
            const std::string &samples = "-",
            const bool &missingToMean = false){

    // initialize stream
    Param param( file, field, region, samples, std::numeric_limits<int>::max(), missingToMean);
    vcfstream vcfObj( param );

    // from VCF, get
    // 1) genotype as arma::mat and 
    // 2) variant properties as varInfo
    // auto [X_geno, vInfo] = vcfObj.getNextChunk();
    DataChunk<arma::mat, VariantInfo> chunk;

    vcfObj.getNextChunk( chunk );
    
    // Convert genotype values for return
    // set colnames as variant IDs
    // set rownames as sample IDs
    VariantInfo info = chunk.getInfo();
    NumericMatrix X = wrap( chunk.getData() );
    colnames(X) = wrap( info.ID );
    rownames(X) = wrap( info.sampleNames );    

    // return genotype data and variant info
    return List::create(    Named("X") = X,
                            Named("info") = toDF(info) );
}


// [[Rcpp::export]]
List extractVcf_eigen( 
            const std::string &file,
            const std::string &field,
            const std::string &region = "",
            const std::string &samples = "-",
            const bool &missingToMean = false){

    // initialize stream

    Param param( file, field, region, samples, std::numeric_limits<int>::max(), missingToMean);
    vcfstream vcfObj( param );

    // from VCF, get
    DataChunk<Eigen::MatrixXd, VariantInfo> chunk;
    vcfObj.getNextChunk( chunk );
    
    // Convert genotype values for return
    // set colnames as variant IDs
    // set rownames as sample IDs
    VariantInfo info = chunk.getInfo();
    NumericMatrix X = wrap( chunk.getData() );
    colnames(X) = wrap( info.ID );
    rownames(X) = wrap( info.sampleNames );    

    // return genotype data and variant info
    return List::create(    Named("X") = X,
                            Named("info") = toDF(info) );
}

// [[Rcpp::export]]
List extractVcf_NM( 
            const std::string &file,
            const std::string &field,
            const std::string &region = "",
            const std::string &samples = "-",
            const bool &missingToMean = false){

    // initialize stream

    Param param( file, field, region, samples, std::numeric_limits<int>::max(), missingToMean);
    vcfstream vcfObj( param );

    // from VCF, get
    DataChunk<Rcpp::NumericMatrix, VariantInfo> chunk;
    vcfObj.getNextChunk( chunk );
    
    // return genotype data and variant info
    return List::create(    Named("X") = chunk.getData(),
                            Named("info") = toDF( chunk.getInfo() ) );
}



// [[Rcpp::export]]
List extractVcf_chunks( 
            const std::string &file,
            const std::string &field,
            const std::string &region = "",
            const std::string &samples = "-",
            const bool &missingToMean = false){

    // initialize stream
    Param param(file, field, region, samples, 2, missingToMean);
    vcfstream vcfObj( param );

    // from VCF, get
    // 1) genotype as arma::mat and 
    // 2) variant properties as varInfo
    // auto [a, b] = vcfObj.getNextChunk();
    // auto [X_geno, vInfo] = vcfObj.getNextChunk();
    DataChunk<arma::mat, VariantInfo> chunk, chunk2;
    vcfObj.getNextChunk( chunk2 );
    vcfObj.getNextChunk( chunk );
    
    // Convert genotype values for return
    // set colnames as variant IDs
    // set rownames as sample IDs
    VariantInfo info = chunk.getInfo();
    NumericMatrix X = wrap( chunk.getData() );
    colnames(X) = wrap( info.ID );
    rownames(X) = wrap( info.sampleNames );    

    // return genotype data and variant info
    return List::create(    Named("X") = X,
                            Named("info") = toDF(info) );
}


// [[Rcpp::export]]
NumericMatrix getDA( const RObject &mat ){

    DelayedStream ds( mat);

    DataChunk<arma::mat, MatrixInfo> chunk;

    ds.getNextChunk( chunk );

    NumericMatrix X = wrap( chunk.getData() );

    return X ;
}



// [[Rcpp::export]]
NumericMatrix getDA_eigen( const RObject &mat ){

    DelayedStream ds( mat);

    DataChunk<Eigen::MatrixXd, MatrixInfo> chunk;

    ds.getNextChunk( chunk );

    NumericMatrix X = wrap( chunk.getData() );

    return X ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix getDA_NM( const RObject &mat ){

    DelayedStream ds( mat);

    DataChunk<Rcpp::NumericMatrix, MatrixInfo> chunk;

    ds.getNextChunk( chunk );

    return chunk.getData();
}








