/***********************************************************************
 * @file        export.cpp
 * @author      Gabriel Hoffman
 * @email       gabriel.hoffman@mssm.edu
 * @brief       Expose GenomicDataStream library to R
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppClock.h>
//[[Rcpp::depends(RcppClock)]


#ifndef DISABLE_EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 


#ifdef _OPENMP
    // [[Rcpp::plugins(openmp)]]
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif


#include <GenomicDataStream.h>

using namespace std;
using namespace vcfpp;
using namespace Rcpp;
using namespace arma;
using namespace gds;

// convert to List
DataFrame toDF( const VariantInfo &vInfo){

    // return created data frame
    return DataFrame::create(
                    Named("CHROM") = Rcpp::wrap(vInfo.CHROM),
                    Named("POS") = Rcpp::wrap(vInfo.POS),
                    Named("ID") = Rcpp::wrap(vInfo.ID),
                    Named("A1") = Rcpp::wrap(vInfo.A1),
                    Named("A2") = Rcpp::wrap(vInfo.A2),
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
    Param param( file, region, samples, std::numeric_limits<int>::max(), missingToMean);
    param.setField(field);
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

    Param param( file, region, samples, std::numeric_limits<int>::max(), missingToMean);
    param.setField(field);
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

    Param param( file, region, samples, std::numeric_limits<int>::max(), missingToMean);    
    param.setField(field);
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

    Param param( file, region, samples, std::numeric_limits<int>::max(), missingToMean);    
    param.setField(field);
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
    Param param(file, region, samples, 2, missingToMean);
    param.setField(field);
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
Rcpp::List getDosage( 
            const std::string &file,
            const std::string &field = "",
            const std::string &region = "",
            const std::string &samples = "-",
            const int &chunkSize = std::numeric_limits<int>::max(),
            const bool &missingToMean = false){

    Param param( file, region, samples, chunkSize, missingToMean);
    param.setField(field);

    unique_ptr<GenomicDataStream> gdsStream = createFileView( param );

    DataChunk<arma::mat, VariantInfo> chunk;

    gdsStream->getNextChunk( chunk );

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




List lm(const Eigen::MatrixXd &X, const Eigen::VectorXd &y) {
    

    Eigen::Index m_p(X.cols());
    Eigen::MatrixXd I_p = Eigen::MatrixXd::Identity(m_p, m_p);
    Eigen::HouseholderQR<Eigen::MatrixXd> QR(X);
    Eigen::VectorXd m_coef      = QR.solve(y);
    // Eigen::VectorXd m_fitted    = X * m_coef;
    Eigen::VectorXd m_se        = QR.matrixQR().topRows(m_p).
        triangularView<Eigen::Upper>().solve(I_p).rowwise().norm();

    return Rcpp::List::create(Rcpp::Named("coefficients") = wrap(m_coef),
                              Rcpp::Named("stderr")       = wrap(m_se),
                              Rcpp::Named("df.residual")  = X.rows() - X.cols());
}

struct ModelFit {
    vec coef;
    vec std_err;
    double df;
    string ID;
    ModelFit() {}
    ModelFit( const vec &coef, const vec &std_err, const double &df) : 
        coef(coef), std_err(std_err), df(df) {}
};


// adapted from https://github.com/RcppCore/RcppArmadillo/blob/master/src/fastLm.cpp
ModelFit lm(const arma::mat& X, const arma::colvec& y) {
    int n = X.n_rows, k = X.n_cols;

    arma::colvec coef = solve(X, y);     // fit model y ~ X
    arma::colvec res  = y - X*coef;            // residuals
    double s2 = dot(res, res) / (n - k); // std.errors of coefficients
    arma::colvec std_err = arma::sqrt(s2 * diagvec(inv(trans(X)*X)));

    return ModelFit( coef, std_err, n-k);
}



vector<ModelFit> linearRegression(const arma::vec &y, const arma::mat &X_cov, const arma::mat &X_features, const VariantInfo &info, const int &nthreads = 1){

    int n_covs = X_cov.n_cols;

    vector<ModelFit> fitList(X_features.n_cols, ModelFit());

    #ifdef _OPENMP 
        // set threads
        omp_set_num_threads(nthreads);
        // disable nested parallelism
        omp_set_max_active_levels(1);
    #endif

    #pragma omp parallel
    {
        // create design matrix with jth feature in the last column
        // X = cbind(X_cov, X_features[,0])
        arma::mat X(X_cov);
        X.insert_cols(n_covs, X_features.col(0));

        // iterate through responses 
        #pragma omp for         
        for(int j=0; j<X_features.n_cols; j++){
            // Create design matrix with intercept as first column
            X.col(n_covs) = X_features.col(j);

            // linear regression        
            ModelFit fit = lm(X, y);
            fit.ID = info.ID[j];

            // save result to list
            fitList.at(j) =  fit;
        }  
    }  

    return fitList;
}

typedef vector<ModelFit> ModelFitList;

List toList( const vector<ModelFitList> & fitList){

    int ncoef = fitList[0][0].coef.n_elem;
    int nrow = 0;
    for(int i=0; i< fitList.size(); i++){
        nrow += fitList[i].size();
    }

    arma::mat coefMat(nrow,ncoef);
    arma::mat seMat(nrow,ncoef);
    arma::vec dfVec(nrow);
    vector<string> ID;
    ID.reserve(nrow);

    int k=0;
    for(int i=0; i< fitList.size(); i++){
        for(int j=0; j< fitList[i].size(); j++){
            coefMat.row(k) = fitList[i][j].coef.t();
            seMat.row(k) = fitList[i][j].std_err.t();
            dfVec(k) = fitList[i][j].df;
            ID.push_back(fitList[i][j].ID);
            k++;
        }
    }

    List lst = List::create(
        Named("ID") = wrap(ID),
        Named("coef") = coefMat,
        Named("se") = seMat,
        Named("df") = dfVec
      );

    return lst;
}



// [[Rcpp::export]]
List fastLM( const arma::colvec& y, 
                const std::string &file,
                const std::string &field = "",
                const std::string &region = "",
                const std::string &samples = "-",
                const int &chunkSize = 4,
                const bool &missingToMean = false, 
                const int &nthreads = 1){

    Param param(file, region, samples, chunkSize, missingToMean);

    param.setField(field);

    // Initialise GenomicDataStream with file
    unique_ptr<GenomicDataStream> gdsStream = createFileView( param );

    if( gdsStream->n_samples() != y.size() ){
        Rcpp::stop("Data stream and y must have same number of samples");
    }

    DataChunk<arma::mat, VariantInfo> chunk;

    // store dosage from chunk
    // n samples and p features
    VariantInfo info_chunk;
    vector<ModelFitList> results;
    arma::mat X_cov(y.n_elem, 1, fill::ones);

    int nVariants = 0;

    while( gdsStream->getNextChunk( chunk ) ){

        // get variant information
        info_chunk = chunk.getInfo();

        // Linear regression with the jth feature
        // used as a covariate in the jth model
        ModelFitList fitList = linearRegression(y, X_cov, chunk.getData(), info_chunk, nthreads);

        nVariants += info_chunk.size();
        // Rcpp::Rcout << "\rVariants processed: " << nVariants << "      ";

        // save results to list
        results.push_back(fitList);
    }
    Rcpp::Rcout << endl;

    return toList( results );
}







/*
// adapted from https://github.com/RcppCore/RcppArmadillo/blob/master/src/fastLm.cpp
List lm(const arma::mat& X, const arma::colvec& y) {
    int n = X.n_rows, k = X.n_cols;

    arma::colvec coef = solve(X, y);     // fit model y ~ X
    arma::colvec res  = y - X*coef;            // residuals
    double s2 = dot(res, res) / (n - k); // std.errors of coefficients
    arma::colvec std_err = arma::sqrt(s2 * diagvec(inv(trans(X)*X)));

    return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                              Rcpp::Named("stderr")       = std_err,
                              Rcpp::Named("df.residual")  = n - k);
}

// adapted from https://github.com/RcppCore/RcppEigen/blob/a0f0564c335316d19b054713642dd5c8bd5084c8/src/fastLm.cpp#L116
List lm(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::VectorXd> &y) {
    
    Eigen::Index m_p(X.cols());
    Eigen::MatrixXd I_p = Eigen::MatrixXd::Identity(m_p, m_p);
    Eigen::HouseholderQR<Eigen::MatrixXd> QR(X);
    Eigen::VectorXd m_coef      = QR.solve(y);
    // Eigen::VectorXd m_fitted    = X * m_coef;
    Eigen::VectorXd m_se        = QR.matrixQR().topRows(m_p).
        triangularView<Eigen::Upper>().solve(I_p).rowwise().norm();

    return Rcpp::List::create(Rcpp::Named("coefficients") = wrap(m_coef),
                              Rcpp::Named("stderr")       = wrap(m_se),
                              Rcpp::Named("df.residual")  = X.rows() - X.cols());
}

List lm(const Eigen::MatrixXd &X, const Eigen::VectorXd &y) {
    

    Eigen::Index m_p(X.cols());
    Eigen::MatrixXd I_p = Eigen::MatrixXd::Identity(m_p, m_p);
    Eigen::HouseholderQR<Eigen::MatrixXd> QR(X);
    Eigen::VectorXd m_coef      = QR.solve(y);
    // Eigen::VectorXd m_fitted    = X * m_coef;
    Eigen::VectorXd m_se        = QR.matrixQR().topRows(m_p).
        triangularView<Eigen::Upper>().solve(I_p).rowwise().norm();

    return Rcpp::List::create(Rcpp::Named("coefficients") = wrap(m_coef),
                              Rcpp::Named("stderr")       = wrap(m_se),
                              Rcpp::Named("df.residual")  = X.rows() - X.cols());
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


List linearRegression(const Eigen::Map<Eigen::VectorXd> &y, const Eigen::Map<Eigen::MatrixXd> &X_cov, const Eigen::Map<Eigen::MatrixXd> &X_features, const VariantInfo &info){

    // create design matrix with jth feature in the last column
    // X = cbind(X_cov, X_features[,0])
    int n_covs = X_cov.cols();
    Eigen::MatrixXd X(X_cov);
    X.conservativeResize(Eigen::NoChange, n_covs+1);
    X.col(n_covs) = X_features.col(0);

    List lst;

    for(int j=0; j<X_features.cols(); j++){
        // Create design matrix with intercept as first column
        X.col(n_covs) = X_features.col(j);

        // linear regression
        List fit = lm(X, y);

        // save result to list
        lst.push_back( fit, info.ID[j] );
    }    

    return lst;
}

List linearRegression(const Eigen::VectorXd &y, const Eigen::MatrixXd &X_cov, const Eigen::MatrixXd &X_features, const VariantInfo &info){

    // create design matrix with jth feature in the last column
    // X = cbind(X_cov, X_features[,0])
    int n_covs = X_cov.cols();
    Eigen::MatrixXd X(X_cov);
    X.conservativeResize(Eigen::NoChange, n_covs+1);
    X.col(n_covs) = X_features.col(0);

    List lst;

    for(int j=0; j<X_features.cols(); j++){
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
List fastLM_eigen( 
                const Eigen::VectorXd& y, 
                const std::string &file,
                const std::string &field,
                const std::string &region = "",
                const std::string &samples = "-",
                const int &chunkSize = 4,
                const bool &missingToMean = false){

    Param param(file, region, samples, chunkSize, missingToMean);

    param.setField(field);

    // Initialise GenomicDataStream with file
    unique_ptr<GenomicDataStream> gdsStream = createFileView( param );

    if( gdsStream->n_samples() != y.size() ){
        Rcpp::stop("Data stream and y must have same number of samples");
    }

    // DataChunk<arma::mat, VariantInfo> chunk;
    DataChunk<Eigen::MatrixXd, VariantInfo> chunk;

    // store dosage from chunk
    // n samples and p features
    VariantInfo info_chunk;
    List lst;
    Eigen::MatrixXd X_cov = Eigen::MatrixXd::Ones(y.size(),1);

    int nVariants = 0;

    while( gdsStream->getNextChunk( chunk ) ){

        // get variant information
        info_chunk = chunk.getInfo();

        // Linear regression with the jth feature
        // used as a covariate in the jth model
        List lst_local = linearRegression(y, X_cov, chunk.getData(), info_chunk);

        nVariants += info_chunk.size();
        Rcpp::Rcout << "nVariants: " << nVariants << endl;

        // save results to list
        append(lst, lst_local);
    }

    return lst;
}


// [[Rcpp::export]]
List fastLM( const arma::colvec& y, 
                const std::string &file,
                const std::string &field,
                const std::string &region = "",
                const std::string &samples = "-",
                const int &chunkSize = 4,
                const bool &missingToMean = false){

    Param param(file, region, samples, chunkSize, missingToMean);

    param.setField(field);

    // Initialise GenomicDataStream with file
    unique_ptr<GenomicDataStream> gdsStream = createFileView( param );

    if( gdsStream->n_samples() != y.size() ){
        Rcpp::stop("Data stream and y must have same number of samples");
    }

    DataChunk<arma::mat, VariantInfo> chunk;

    // store dosage from chunk
    // n samples and p features
    VariantInfo info_chunk;
    List lst, lst_store;
    arma::mat X_cov(y.n_elem, 1, fill::ones);

    int nVariants = 0;

    Rcpp::Clock clock;

    while( gdsStream->getNextChunk( chunk ) ){

        // get variant information
        info_chunk = chunk.getInfo();

        // Linear regression with the jth feature
        // used as a covariate in the jth model

        clock.tick("linearRegression");
        List lst_local = linearRegression(y, X_cov, chunk.getData(), info_chunk);
        clock.tock("linearRegression");

        nVariants += info_chunk.size();
        Rcpp::Rcout << "nVariants: " << nVariants << endl;

        // save results to list
        clock.tick("append");
        // append(lst, lst_local);
        lst_store.push_back(lst_local);
        clock.tock("append");
    }

    clock.stop("info");

    return lst_store;
}


*/


