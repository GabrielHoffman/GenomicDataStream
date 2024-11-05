/***********************************************************************
 * @file		regerssion.cpp
 * @author	   	Gabriel Hoffman
 * @email	   	gabriel.hoffman@mssm.edu
 * @brief	   	Expose regression functions to R
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


// DISABLE warning: solve(): system is singular
#define ARMA_WARN_LEVEL 1

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
 
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


struct ModelFit {
	vec coef;
	vec se;
	double df;
	string ID;
    ModelFit() {}
	ModelFit( const vec &coef, const vec &se, const double &df) : 
		coef(coef), se(se), df(df) {}
};



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
            seMat.row(k) = fitList[i][j].se.t();
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




/** Linear regression 
 * 
 * @param X design matrix
 * @param y response vector
 * 
/ adapted from https://github.com/RcppCore/RcppArmadillo/blob/master/src/fastLm.cpp
/ https://genomicsclass.github.io/book/pages/qr_and_regression.html
*/
ModelFit lm(const arma::mat& X, const arma::colvec& y) {

	int n = X.n_rows, k = X.n_cols;

    vec beta = solve(X, y);     // fit model y ~ X
    vec res  = y - X*beta;            // residuals
    double s2 = dot(res, res) / (n - k); // std.errors of coefficients
    vec se = sqrt(s2 * diagvec(pinv(trans(X)*X)));

	return ModelFit( beta, se, n-k);
}

// [[Rcpp::export]]
List test_lm(const arma::mat& X, const arma::colvec& y) {

	ModelFit fit = lm(X,y);

	return Rcpp::List::create(Rcpp::Named("coefficients")   = fit.coef,
							  Rcpp::Named("se")			 = fit.se);
}	

vector<ModelFit> linearRegression(const arma::vec &y, const arma::mat &X_cov, const arma::mat &X_features, const DataInfo *info, const int &nthreads = 1){

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
			fit.ID = info->getFeatureName(j);

			// save result to list
			fitList.at(j) =  fit;
		}  
	}  

	return fitList;
}


vector<ModelFit> linearRegressionResponses(const arma::mat &Y, const arma::mat &X, const DataInfo *info, const int &nthreads = 1){

    int n_covs = X.n_cols;

    vector<ModelFit> fitList(Y.n_cols, ModelFit());

    #ifdef _OPENMP 
        // set threads
        omp_set_num_threads(nthreads);
        // disable nested parallelism
        omp_set_max_active_levels(1);
    #endif

    #pragma omp parallel
    {
        // iterate through responses 
        #pragma omp for      
        for(int j=0; j<Y.n_cols; j++){

            // linear regression        
            ModelFit fit = lm(X, Y.col(j));
            fit.ID = info->getFeatureName(j);

            // save result to list
            fitList.at(j) =  fit;
        }  
    }  

    return fitList;
}

// [[Rcpp::export]]
List fastLM( const arma::colvec& y, 
				const std::string &file,
				const std::string &field = "",
				const std::string &region = "",
				const std::string &samples = "-",
				const int &chunkSize = 4,
				const bool &missingToMean = false, 
				const int &nthreads = 1,
                const bool &verbose = true){

	Param param(file, region, samples, chunkSize, missingToMean);
	param.setField(field);

	// Initialise GenomicDataStream with file
	unique_ptr<GenomicDataStream> gdsStream = createFileView( param );

	if( gdsStream->n_samples() != y.size() ){
		Rcpp::stop("Data stream and y must have same number of samples");
	}

	DataChunk<arma::mat> chunk;

	// store dosage from chunk
	// n samples and p features
	VariantInfo *info_chunk;
	vector<ModelFitList> results;
	arma::mat X_cov(y.n_elem, 1, fill::ones);

	int nModels = 0;

	while( gdsStream->getNextChunk( chunk ) ){

		// get variant information
		info_chunk = chunk.getInfo<VariantInfo>();

		// Linear regression with the jth feature
		// used as a covariate in the jth model
		ModelFitList fitList = linearRegression(y, X_cov, chunk.getData(), info_chunk, nthreads);

		nModels += info_chunk->size();

		if( verbose ) 
            Rcpp::Rcout << "\rModels fit: " << nModels << "	  ";

		// save results to list
		results.push_back(fitList);
	}
	if( verbose ) Rcpp::Rcout << endl;

	return toList( results );
}


// [[Rcpp::export]]
List regrExprResponse(
                const RObject &mat, 
                const vector<string> &rowNames, 
                const int &chunkSize,
                const int &nthreads,
                const bool &verbose = true ){

    DelayedStream ds( mat, rowNames, chunkSize);

    DataChunk<arma::mat> chunk;
    MatrixInfo *info;
    vector<ModelFitList> results;
    arma::mat X_design( ds.n_samples(), 1, fill::ones);

    int nModels = 0;

    while( ds.getNextChunk( chunk ) ){

        // get variant information
        info = chunk.getInfo<MatrixInfo>();

        // Linear regression with the jth feature
        // used as a covariate in the jth model
        ModelFitList fitList = linearRegressionResponses(chunk.getData(), X_design, info, nthreads);

        nModels += info->size();

        if( verbose ) 
            Rcpp::Rcout << "\rModels fit: " << nModels << "   ";

        // save results to list
        results.push_back(fitList);
    }
    if( verbose ) Rcpp::Rcout << endl;

    return toList( results );
}


