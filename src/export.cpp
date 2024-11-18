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


#include <GenomicDataStream.h>

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
	DataChunk<arma::mat> chunk;

	vcfObj.getNextChunk( chunk );
	
	// Convert genotype values for return
	// set colnames as variant IDs
	// set rownames as sample IDs
	VariantInfo *info = chunk.getInfo<VariantInfo>();
	NumericMatrix X = wrap( chunk.getData() );
	colnames(X) = wrap( info->getFeatureNames() );
	rownames(X) = wrap( info->sampleNames );	

	// return genotype data and variant info
	return List::create(	Named("X") = X,
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
	DataChunk<Eigen::MatrixXd> chunk;
	vcfObj.getNextChunk( chunk );
	
	// Convert genotype values for return
	// set colnames as variant IDs
	// set rownames as sample IDs
    VariantInfo *info = chunk.getInfo<VariantInfo>();
	NumericMatrix X = wrap( chunk.getData() );
	colnames(X) = wrap( info->getFeatureNames() );
	rownames(X) = wrap( info->sampleNames );	

	// return genotype data and variant info
	return List::create(	Named("X") = X,
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
	DataChunk<Rcpp::NumericMatrix> chunk;
	vcfObj.getNextChunk( chunk );
	
	// return genotype data and variant info
	return List::create(	Named("X") = chunk.getData(),
							Named("info") = toDF( chunk.getInfo<VariantInfo>() ) );
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
	DataChunk<vector<double>> chunk;
	vcfObj.getNextChunk( chunk );
	
	// return genotype data and variant info
	return List::create(	Named("X") = Rcpp::wrap(chunk.getData()),
							Named("info") = toDF( chunk.getInfo<VariantInfo>() ) );
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
	DataChunk<arma::mat> chunk, chunk2;
	vcfObj.getNextChunk( chunk2 ); 
	vcfObj.getNextChunk( chunk );
	
	// Convert genotype values for return
	// set colnames as variant IDs
	// set rownames as sample IDs
	VariantInfo *info = chunk.getInfo<VariantInfo>();
	NumericMatrix X = wrap( chunk.getData() );
	colnames(X) = wrap( info->getFeatureNames() );
	rownames(X) = wrap( info->sampleNames );	

	// return genotype data and variant info
	return List::create(	Named("X") = X,
							Named("info") = toDF(info) );
}


CharacterVector convert_to_character_vector(const vector<string>& vec) {
  CharacterVector result(vec.size());

  for (size_t i = 0; i < vec.size(); ++i) {
    result[i] = vec[i];
  }

  return result;
}

// [[Rcpp::export(rng=false)]]
List getDA( const RObject &mat, const vector<string> &rowNames, const int &chunkSize ){

    DelayedStream ds( mat, rowNames, chunkSize);

    DataChunk<arma::mat> chunk;

    ds.getNextChunk( chunk );
    // ds.getNextChunk( chunk );

    MatrixInfo *info = chunk.getInfo<MatrixInfo>();

    NumericMatrix X = wrap(chunk.getData());

    // return genotype data and variant info
    return List::create(    Named("X") = X,
                            Named("info") = convert_to_character_vector(info->getFeatureNames()) );
}









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

	DataChunk<arma::mat> chunk;

	gdsStream->getNextChunk( chunk );

	// Convert genotype values for return
	// set colnames as variant IDs
	// set rownames as sample IDs
	VariantInfo *info = chunk.getInfo<VariantInfo>();
	NumericMatrix X = wrap( chunk.getData() );
	colnames(X) = wrap( info->getFeatureNames() );
	rownames(X) = wrap( info->sampleNames );	

	// return genotype data and variant info
	return List::create(	Named("X") = X,
							Named("info") = toDF(info) );
}



// [[Rcpp::export]]
SEXP create_xptr( 
			const std::string &file,
			const std::string &field = "",
			const std::string &region = "",
			const std::string &samples = "-",
			const int &chunkSize = std::numeric_limits<int>::max(),
			const bool &missingToMean = false){

	Param param( file, region, samples, chunkSize, missingToMean);
	param.setField(field);
 
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
	// returns true is chunk is valid
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




