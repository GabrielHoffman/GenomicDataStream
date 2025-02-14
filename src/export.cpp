/***********************************************************************
 * @file		export.cpp
 * @author	   	Gabriel Hoffman
 * @email	   	gabriel.hoffman@mssm.edu
 * @brief	   	Expose GenomicDataStream library to R
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
#include "DataTable.h"

using namespace std;
using namespace vcfpp;
using namespace Rcpp;
using namespace arma;
using namespace gds;

// convert to List
DataFrame toDF( const VariantInfo *vInfo){

	// return created data frame
	return DataFrame::create(
					Named("CHROM") = Rcpp::wrap(vInfo->CHROM),
					Named("POS") = Rcpp::wrap(vInfo->POS),
					Named("ID") = Rcpp::wrap(vInfo->getFeatureNames()),
					Named("A1") = Rcpp::wrap(vInfo->A1),
					Named("A2") = Rcpp::wrap(vInfo->A2),
					_["stringsAsFactors"] = false);
}



// [[Rcpp::export]]
SEXP create_xptr( 
			const std::string &file,
			const std::string &field = "",
			const std::string &region = "",
			const std::string &samples = "-",
			const double &minVariance = 0,
			const int &chunkSize = 10000,
			const bool &missingToMean = true){

	Param param( file, region, samples, minVariance, chunkSize, missingToMean);
	param.setField(field);
 
 	// calls constructor for GenomicDataStream
	Rcpp::XPtr<BoundDataStream> z( new BoundDataStream(param), true);

	return z;
}


// [[Rcpp::export]]
List getInfo(SEXP x){

	Rcpp::XPtr<BoundDataStream> ptr(x);

	return List::create(Named("streamType") = ptr->ptr->getStreamType(),
						Named("nsamples") = ptr->ptr->n_samples());
}


// [[Rcpp::export]]
SEXP setRegions_rcpp( SEXP x, const string &regionString){
	Rcpp::XPtr<BoundDataStream> ptr(x);

	vector<string> regions = splitRegionString( regionString );

	ptr->ptr->setRegions( regions );

	ptr->atEndOfStream = false;
	ptr->featuresRead = 0;

	return ptr;
}


// [[Rcpp::export]]
CharacterVector getSampleNames_rcpp( SEXP x){
	Rcpp::XPtr<BoundDataStream> ptr(x);

	vector<string> IDs = ptr->ptr->getSampleNames();

	return wrap(IDs);
}


// [[Rcpp::export]]
bool atEndOfStream_rcpp( SEXP x){
	Rcpp::XPtr<BoundDataStream> ptr(x);
	return ptr->atEndOfStream;
}

// [[Rcpp::export]]
long featuresRead_rcpp( SEXP x){
	Rcpp::XPtr<BoundDataStream> ptr(x);
	return ptr->featuresRead;
}


// [[Rcpp::export]]
List getNextChunk_rcpp( SEXP x){ 

	Rcpp::XPtr<BoundDataStream> ptr(x);

	DataChunk<arma::mat> chunk;

	// get chunk of data, 
	// returns true if chunk is valid
	bool isValid = ptr->ptr->getNextChunk( chunk );

	// if not valid
	// set atEndOfStream to true and return empty list
	if( ! isValid ){
		ptr->atEndOfStream = true;
		return List::create();
	}

	// else continue with valid data

	// Convert genotype values for return
	// set colnames as variant IDs
	// set rownames as sample IDs
	VariantInfo *info = chunk.getInfo<VariantInfo>();
	NumericMatrix X = wrap( chunk.getData() );
	colnames(X) = wrap( info->getFeatureNames() );
	rownames(X) = wrap( info->sampleNames );	

	ptr->featuresRead += info->size();

	// return genotype data and variant info
	return List::create(	Named("X") = X,
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

// [[Rcpp::export]]
void test_DataTable(const string &file, const string &headerKey, const char delim='\t'){

	DataTable dt(file, headerKey, delim);

	dt.print(Rcpp::Rcout, "\t");
}


