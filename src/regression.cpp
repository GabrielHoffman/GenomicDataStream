/***********************************************************************
 * @file		regression.cpp
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


#include "GenomicDataStream.h"
#include <linearRegression.h>
#include <glm.h>
#include <exportToR.h>

using namespace std;
using namespace vcfpp;
using namespace Rcpp;
using namespace arma;
using namespace gds;

// Concatenate then run toList() on combined vector
const List toList( const vector<ModelFitList> & fitList){

	// count number of models
	int nmodels = 0;
    for(int k=0; k < fitList.size(); k++){
        nmodels += fitList[k].size();
    }

    // allocate memory
	vector<ModelFit> fitListCombined;
	fitListCombined.reserve( nmodels );

 	for(int k=0; k < fitList.size(); k++){
        fitListCombined.insert( fitListCombined.end(), 
        						fitList[k].begin(),
        						fitList[k].end());
    }

    return toList( fitListCombined );
}


// Concatenate then run toList() on combined vector
const List toList( const vector<GLMModelFitList> & fitList){

	// count number of models
	int nmodels = 0;
    for(int k=0; k < fitList.size(); k++){
        nmodels += fitList[k].size();
    }

    // allocate memory
	vector<GLMModelFit> fitListCombined;
	fitListCombined.reserve( nmodels );

 	for(int k=0; k < fitList.size(); k++){
        fitListCombined.insert( fitListCombined.end(), 
        						fitList[k].begin(),
        						fitList[k].end());
    }

    return toList( fitListCombined );
}



// [[Rcpp::export]]
List test_lm(const arma::mat& X, const arma::colvec& y) {

	ModelFit fit = lm(X,y);

	return Rcpp::List::create(Rcpp::Named("coefficients")   = fit.coef,
							  Rcpp::Named("se")			 = fit.se);
}	




// [[Rcpp::export]]
List lmFitFeatures_export( const arma::colvec& y, 
					const arma::mat &X_design,
					List gds, 					
					const arma::vec &weights,
					const int &detail = 0,
					const bool &preprojection=true, 
					const int &nthreads=1){

	// Get parameters from R::GenomicDataStream
	string file = gds["file"];
	string field = gds["field"];
	string region = gds["region"];
	string samples = gds["samples"];
	int chunkSize = gds["chunkSize"];
	bool missingToMean = gds["missingToMean"];

	// Initialize C++-level GenomicDataStream
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

	int nModels = 0;
	ModelDetail md = LOW;

	while( gdsStream->getNextChunk( chunk ) ){

		// get variant information
		info_chunk = chunk.getInfo<VariantInfo>();

		// Linear regression with the jth feature
		// used as a covariate in the jth model
		ModelFitList fitList = lmFitFeatures(y, X_design, chunk.getData(), info_chunk->getFeatureNames(), weights, md, preprojection, nthreads);

		nModels += fitList.size();

		// save results to list
		results.push_back(fitList);
	}

	return toList( results );
}




// [[Rcpp::export]]
List lmFitResponses_export(
                const RObject &mat, 
                const arma::mat & X_design,
                const vector<string> &ids,
				const arma::mat &Weights,
                const int &chunkSize,
				const int &detail = 0, 
                const int &nthreads = 1,
                const bool &verbose = true ){

    DelayedStream ds( mat, ids, chunkSize);

    DataChunk<arma::mat> chunk;
    MatrixInfo *info;
    vector<ModelFitList> results;

    ModelDetail md = LOW;
    int nModels = 0;
    arma::mat W;

    while( ds.getNextChunk( chunk ) ){

        // get variant information
        info = chunk.getInfo<MatrixInfo>();

        // weights for this chunk
        W = Weights.rows(nModels, nModels+info->size()-1).t();

        // Linear regression with the jth feature
        // used as a covariate in the jth model
        ModelFitList fitList = lmFitResponses(chunk.getData(), X_design, info->getFeatureNames(), W, md, nthreads);

        nModels += fitList.size();

        if( verbose ) 
            Rcpp::Rcout << "\rModels fit: " << nModels << "   ";

        // save results to list
        results.push_back(fitList);
    }
    if( verbose ) Rcpp::Rcout << endl;

    return toList( results );
}


// [[Rcpp::export]]
List glmFitFeatures_export( const arma::colvec& y, 
							const arma::mat &X_design,
							List gds,
							const std::string &family, 
							const arma::vec &weights, 
							const arma::vec &offset,
							const int detail = 1, 
							const bool &doCoxReid = true, 
							const bool &shareTheta = false,
							const bool &fastApprox = false,
							const int &nthreads = 1, 
							const double &epsilon = 1e-8, 
							const double &maxit = 25, 
							const double &epsilon_nb = 1e-4, 
							const double & maxit_nb = 5	){

	// Get parameters from R::GenomicDataStream
	string file = gds["file"];
	string field = gds["field"];
	string region = gds["region"];
	string samples = gds["samples"];
	int chunkSize = gds["chunkSize"];
	bool missingToMean = gds["missingToMean"];

	// Initialize C++-level GenomicDataStream
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
	vector<GLMModelFitList> results;

	int nModels = 0;
	ModelDetail md = LOW;

	while( gdsStream->getNextChunk( chunk ) ){

		// get variant information
		info_chunk = chunk.getInfo<VariantInfo>();

		// Linear regression with the jth feature
		// used as a covariate in the jth model
		GLMModelFitList fitList = glmFitFeatures(y, X_design, chunk.getData(), info_chunk->getFeatureNames(), family, weights, offset, md, doCoxReid, shareTheta, fastApprox, nthreads, epsilon, maxit, epsilon_nb, maxit_nb);

		nModels += fitList.size();

		// save results to list
		results.push_back(fitList);
	}

	return toList( results );
}

// [[Rcpp::export]]
List glmFitResponses_export(
                const RObject &mat, 
                const arma::mat & X_design,
                const vector<string> &ids,
				const vector<string> &family, 
				const arma::vec &weights, 
				const arma::vec &offset,
                const int &chunkSize,
				const int detail = 1, 
				const bool &doCoxReid = true, 
				const int &nthreads = 1, 
                const double &epsilon = 1e-8, 
                const double &maxit = 25, 
                const double &epsilon_nb = 1e-4, 
                const double & maxit_nb = 5){

    DelayedStream ds( mat, ids, chunkSize);

    DataChunk<arma::mat> chunk;
    MatrixInfo *info;
    vector<GLMModelFitList> results;

    ModelDetail md = LOW;
    int nModels = 0;

    while( ds.getNextChunk( chunk ) ){

        // get variant information
        info = chunk.getInfo<MatrixInfo>();

        // family values for this chunk
        auto start = family.begin() + nModels;
        auto end = start + info->size();
        vector<string> fam(start, end);

        // Linear regression with the jth feature
        // used as a covariate in the jth model
        GLMModelFitList fitList = glmFitResponses(chunk.getData(), X_design, info->getFeatureNames(), fam, weights, offset, md, doCoxReid, nthreads, epsilon, maxit, epsilon_nb, maxit_nb);

        nModels += fitList.size();

        // save results to list
        results.push_back(fitList);
    }

    return toList( results );
}

