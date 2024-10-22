/***********************************************************************
 * @file        export.cpp
 * @author      Gabriel Hoffman
 * @email       gabriel.hoffman@mssm.edu
 * @brief       Expose GenomicDataStream library to R
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#ifdef USE_EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 

#include <GenomicDataStream.h>

using namespace std;
using namespace vcfpp;
using namespace Rcpp;
using namespace arma;
using namespace gds;


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
List extractVcf_vector( 
            const std::string &file,
            const std::string &field,
            const std::string &region = "",
            const std::string &samples = "-",
            const bool &missingToMean = false){

    // initialize stream

    Param param( file, field, region, samples, std::numeric_limits<int>::max(), missingToMean);
    vcfstream vcfObj( param );

    // from VCF, get
    DataChunk<vector<double>, VariantInfo> chunk;
    vcfObj.getNextChunk( chunk );
    
    // return genotype data and variant info
    return List::create(    Named("X") = Rcpp::wrap(chunk.getData()),
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


// // [[Rcpp::export]]
// NumericMatrix getDA( const RObject &mat ){

//     DelayedStream ds( mat);

//     DataChunk<arma::mat, MatrixInfo> chunk;

//     ds.getNextChunk( chunk );

//     NumericMatrix X = wrap( chunk.getData() );

//     return X ;
// }



// // [[Rcpp::export]]
// NumericMatrix getDA_eigen( const RObject &mat ){

//     DelayedStream ds( mat);

//     DataChunk<Eigen::MatrixXd, MatrixInfo> chunk;

//     ds.getNextChunk( chunk );

//     NumericMatrix X = wrap( chunk.getData() );

//     return X ;
// }

// // [[Rcpp::export]]
// Rcpp::NumericMatrix getDA_NM(  RObject mat ){


//     Rcpp::Rcout << "read_lin_block" << std::endl;
//     auto a = beachmat::read_lin_block(mat);
//     Rcpp::Rcout << "end" << std::endl;
    
//     Rcpp::Rcout << "DelayedStream" << std::endl;
//     DelayedStream ds( mat);
//     Rcpp::Rcout << "success" << std::endl;

//     DataChunk<Rcpp::NumericMatrix, MatrixInfo> chunk;


//     Rcpp::Rcout << "getNextChunk" << std::endl;
//     ds.getNextChunk( chunk );
//     Rcpp::Rcout << "success" << std::endl;

//     return chunk.getData();
// }


// // [[Rcpp::export]]
// Rcpp::NumericVector getDA_vector(const RObject &mat ){

//     DelayedStream ds( mat);

//     DataChunk<vector<double>, MatrixInfo> chunk;

//     ds.getNextChunk( chunk );

//     return Rcpp::wrap(chunk.getData());
// }


// #include "Rtatami.h"
// #include <algorithm>

// // [[Rcpp::export]]
// Rcpp::NumericVector column_sums(const Rcpp::RObject &initmat) {
//     Rtatami::BoundNumericPointer parsed(initmat);
//     const auto& ptr = parsed->ptr;

//     auto NR = ptr->nrow();
//     auto NC = ptr->ncol();
//     std::vector<double> buffer(NR);
//     Rcpp::NumericVector output(NC);
//     auto wrk = ptr->dense_column();

//     for (int i = 0; i < NC; ++i) {
//         auto extracted = wrk->fetch(i, buffer.data());
//         output[i] = std::accumulate(extracted, extracted + NR, 0.0);
//     }

//     DelayedStream ds(initmat);


//     return output;
// }



// [[Rcpp::export]]
Rcpp::List test_bgen( 
            const std::string &file,
            const std::string &field,
            const std::string &region = "",
            const std::string &samples = "-",
            const int &chunkSize = std::numeric_limits<int>::max(),
            const bool &missingToMean = false){

    Param param( file, field, region, samples, chunkSize, missingToMean);

    bgenstream bgenObj(param);

    DataChunk<arma::mat, VariantInfo> chunk;

    bgenObj.getNextChunk( chunk );

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
arma::vec colSums_test( const arma::mat &X){
    return colSums(X);
}

// [[Rcpp::export]]
void standardize_test( arma::mat &X, const bool &center = true, const bool &scale = true ){

    standardize(X, center, scale);
}


// adapted from https://github.com/RcppCore/RcppArmadillo/blob/master/src/fastLm.cpp
List lm(const arma::mat& X, const arma::colvec& y) {
    int n = X.n_rows, k = X.n_cols;

    arma::colvec coef = arma::solve(X, y);     // fit model y ~ X
    arma::colvec res  = y - X*coef;            // residuals
    double s2 = arma::dot(res, res) / (n - k); // std.errors of coefficients
    arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));

    return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                              Rcpp::Named("stderr")       = std_err,
                              Rcpp::Named("df.residual")  = n - k);
}

List linearRegression(const arma::vec &y, const arma::mat &X_cov, const arma::mat &X_features, const VariantInfo &info){

    // create design matrix with jth feature in the last column
    // X = cbind(X_cov, X_features[,0])
    int n_covs = X_cov.n_cols;
    arma::mat X(X_cov);
    X.insert_cols(n_covs, X_features.col(0));

    List lst;

    for(int j=0; j<X_features.n_cols; j++){
        // Create design matrix with intercept as first column
        X.col(n_covs) = X_features.col(j);

        // linear regression
        List fit = lm(X, y);

        // save result to list
        lst.push_back( fit, info.ID[j] );
    }    

    return lst;
}


void append(List &lst, const List &lst2){
    CharacterVector ch = lst2.names();
    // append entries to List
    for(int i=0; i<lst2.size(); i++){
        lst.push_back(lst2[i], as<string>(ch[i]));
    }
}


// [[Rcpp::export]]
List fastLM( const arma::colvec& y, 
                const std::string &file,
                const std::string &field,
                const std::string &region = "",
                const std::string &samples = "-",
                const bool &missingToMean = false){

    int chunkSize = 4;
    Param param(file, field, region, samples, chunkSize, missingToMean);

    // Initialise GenomicDataStream with file
    unique_ptr<GenomicDataStream> gdsStream = createFileView( param );

    if( gdsStream->n_samples() != y.size() ){
        Rcpp::stop("Data stream and y must have same number of samples");
    }

    DataChunk<arma::mat, VariantInfo> chunk;

    // store dosage from chunk
    // n samples and p features
    arma::mat X_chunk;
    VariantInfo info_chunk;
    List lst;
    arma::mat X_cov(60, 1, fill::ones);

    while( gdsStream->getNextChunk( chunk ) ){

        // get data from chunk
        X_chunk = chunk.getData();

        // get variant information
        info_chunk = chunk.getInfo();

        // standardize X for mean and sd
        // standardize( X_chunk );

        // Linear regression with the jth feature
        // used as a covariate in the jth model
        List lst_local = linearRegression(y, X_cov, X_chunk, info_chunk);

        // save results to list
        append(lst, lst_local);
    }

    return lst;
}


// [[Rcpp::export]]
List test_bgen2( const arma::colvec y, 
                const std::string &file,
                const std::string &field,
                const std::string &region = "",
                const std::string &samples = "-",
                const int &chunkSize = std::numeric_limits<int>::max(),
                const bool &missingToMean = true){ 

    Param param(file, field, region, samples, chunkSize, missingToMean);

    bgenstream bgenObj(param);

    if( bgenObj.n_samples() != y.size() ){
        Rcpp::stop("Data stream and y must have same number of samples");
    }

    DataChunk<arma::mat, VariantInfo> chunk;

    // store dosage from chunk
    // n samples and p features
    arma::mat X_chunk;
    VariantInfo info_chunk;
    List lst;
    arma::mat X_cov(500, 1, fill::ones);

    int i = 0;
    while( bgenObj.getNextChunk( chunk ) ){

        // get data from chunk
        X_chunk = chunk.getData();

        // get variant information
        info_chunk = chunk.getInfo();

        // standardize X for mean and sd
        // standardize( X_chunk );

        // Linear regression with the jth feature
        // used as a covariate in the jth model
        List lst_local = linearRegression(y, X_cov, X_chunk, info_chunk);

        // // save results to list
        append(lst, lst_local);
    }

    return lst;
}

