/***********************************************************************
 * @file		export.cpp
 * @author	  Gabriel Hoffman
 * @email	   gabriel.hoffman@mssm.edu
 * @brief	   Expose GenomicDataStream library to R
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
typedef struct BoundDataStream {
	BoundDataStream(const Param &param){
		ptr = createFileView_shared( param );
	}

	shared_ptr<gds::GenomicDataStream> ptr;
	bool reachedEnd = false;
	long featuresRead = 0;
} BoundDataStream;

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
bool hasReachedEnd_rcpp( SEXP x){
	Rcpp::XPtr<BoundDataStream> ptr(x);
	return ptr->reachedEnd;
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

	ptr->reachedEnd = ! ptr->ptr->getNextChunk( chunk );

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




