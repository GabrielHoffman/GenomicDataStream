/***********************************************************************
 * @file		export.cpp
 * @author	Gabriel Hoffman
 * @email	  gabriel.hoffman@mssm.edu
 * @brief	   Expose GenomicDataStream library to R
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


// DISABLE warning: solve(): system is singular
#include <string>
#define ARMA_WARN_LEVEL 1

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
 
#ifndef DISABLE_EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 

#include "GenomicDataStream.h"
#include "DataTable.h"
#include "ParallelGenomicChunks.h"
#include "PCAOne.h"

#include "Rand.hpp"

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
      const double &MAF = 0,
			const double &minVariance = 0,
			const int &chunkSize = 10000,
			const bool &missingToMean = true){

	Param param( file, region, samples, MAF, minVariance, chunkSize, missingToMean);
	param.setField(field);
 
 	// calls constructor for GenomicDataStream
	Rcpp::XPtr<BoundDataStream> z( new BoundDataStream(param), true);

	return z;
}


// [[Rcpp::export]]
List getInfo(SEXP x){

	Rcpp::XPtr<BoundDataStream> ptr(x);

	return List::create(
    Named("streamType") = ptr->ptr->getStreamType(),
		Named("nsamples")   = ptr->ptr->n_samples());
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
SEXP setChunkSize_rcpp( SEXP x, const double &chunkSize){
	Rcpp::XPtr<BoundDataStream> ptr(x);

	ptr->ptr->setChunkSize( chunkSize );

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

  DataFrame df = toDF(info);

	// return genotype data and variant info
	return List::create(	Named("X") = X,
							Named("info") = df );
}


List summarizeChunks( const shared_ptr<GenomicDataStream> gds){

  DataChunk<arma::mat> chunk;
  VariantInfo *info;
  vector<string> variantIDs, intervals, tmp, regions;
  vector<int> chunkCounts;

  // loop through chunks
  while( gds->getNextChunk( chunk ) ){

      // get variant information
      info = chunk.getInfo<VariantInfo>();

      tmp = info->getFeatureNames();
      variantIDs.insert(variantIDs.end(), tmp.begin(), tmp.end());

      tmp = info->getRegions();
      regions.insert(regions.end(), tmp.begin(), tmp.end());

      intervals.push_back( info->getInterval() );
      chunkCounts.push_back( info->size() );
  }

  // MultiGenomicStreamProcessor<Eigen::MatrixXd> processor(gds->getParam(), regions, nthreads, nthreads);
  // std::mutex pcaMutex;  

  // processor.processChunk([&](const gds::DataChunk<Eigen::MatrixXd> &chunk, size_t b) {
      
  //   std::lock_guard<std::mutex> lock(pcaMutex);

  //   // get variant information
  //   VariantInfo *info = chunk.getInfo<VariantInfo>(); 

  //   tmp = info->getFeatureNames();
  //   variantIDs.insert(variantIDs.end(), tmp.begin(), tmp.end());

  //   tmp = info->getRegions();
  //   regions.insert(regions.end(), tmp.begin(), tmp.end());

  //   intervals.push_back( info->getInterval() );
  //   chunkCounts.push_back( info->size() ); 
  // });


  return List::create(
        Named("intervals") = intervals,
        Named("chunks") = chunkCounts,
        Named("sampleIDs") = gds->getSampleNames(),
        Named("variantIDs") = variantIDs,
        Named("regions") = regions);
}


// [[Rcpp::export]]
List summarizeChunks_rcpp( SEXP x){ 

  Rcpp::XPtr<BoundDataStream> ptr(x);

  return summarizeChunks( ptr->ptr );
}



// [[Rcpp::export]]
arma::vec colSums_test( const arma::mat &X){
	return colSums(X);
}



// [[Rcpp::export(.standardize_in_place)]]
void standardize_in_place( arma::mat &X, const bool &center = true, const bool &scale = true ){

	if( X.n_rows != 0 && X.n_cols != 0){
		standardize(X, center, scale);
	}
}

// [[Rcpp::export]]
void test_DataTable(const string &file, const string &headerKey, const string &delim="\t "){

	DataTable dt(file, headerKey, delim);

	dt.print(Rcpp::Rcout, "\t");
}


// [[Rcpp::export]]
List stream_pcaone( SEXP x, 
                    const string &region,
                    int m, 
                    int k, 
                    int nchunks,
                    int s = 20, 
                    int p = 7, 
                    int B = 64, 
                    int threads = 4,
                    int threads_eigen = 1,
                    const bool verbose = true,
                    const bool scaleAndCenter = true) {

  XPtr<BoundDataStream> ptr(x);
  shared_ptr<GenomicDataStream> gds = ptr->ptr;

  Eigen::setNbThreads(threads_eigen);   

  PCA res = pcaone( gds, region, m, k, nchunks, s, p, B, threads, verbose, scaleAndCenter);

  return List::create(
        Named("d") = wrap(res.d),
        Named("u") = wrap(res.u),
        Named("v") = wrap(res.v),
        Named("featureIds") = wrap(res.featureIds));
}



// [[Rcpp::export]]
Rcpp::List stream_pcaone_robj(
										const RObject &x, 
										const std::vector<std::string> &ids, 
    								const int &n, 
    								const int &chunkSize,
    								const int &nchunks,
										int m, 
										int k, 
										int s = 20, 
										int p = 7, 
										int B = 64, 
										int threads = 4, 
                    int threads_eigen = 1,
										const bool verbose = true,
                    const bool scaleAndCenter = true) {

  // tatami interface to matrix x 
  Rtatami::BoundNumericPointer *parsed = new Rtatami::BoundNumericPointer(x);
  const auto& ptr = (*parsed)->ptr;

	const int l = min(m, k + s);
  size_t band = ceil((double) nchunks / B);
  std::mutex pcaMutex;  

  // initialize global matricies
  auto randomEngine = std::default_random_engine{};
  Eigen::MatrixXd Omg = StandardNormalRandom<Eigen::MatrixXd, std::default_random_engine>(n, l, randomEngine);
  Eigen::MatrixXd Omg2 = Omg;
  Eigen::MatrixXd H1 = Eigen::MatrixXd::Zero(n, l);
  Eigen::MatrixXd H2 = Eigen::MatrixXd::Zero(n, l);
  Eigen::MatrixXd H(n, l), G(m, l), R(l, l), Rt(l, l);

  Eigen::setNbThreads(threads_eigen);

  // power iterations
  for (int pi = 0; pi <= p; pi++) {

    if( verbose ){
      Rcpp::Rcout << "\rEpoch " << pi << " / " << p << "  ";
    }

    if (std::pow(2, pi) >= B) {
      // reset H1, H2 to zero
      H1.setZero();
      H2.setZero();
    }
    band = std::fmin(band * 2, nchunks);

    size_t i{1}, b{0};
    // parallelize across threads
    tatami::parallelize([&](size_t thread_id, int chk, int len) -> void {

      // initialize variables for this thread
      DelayedStream ds( ptr, ids, chunkSize);
      DataChunk<Eigen::MatrixXd> chunk;
      Eigen::MatrixXd Ab;

      // In this thread, loop through multiple chunks
      for(int idx=chk; idx<chk+len; idx++){

        int start = idx*chunkSize;
        ds.getNextChunk( chunk, start, chunkSize);

        Ab = chunk.getData();

        if( scaleAndCenter ){
          standardize(Ab);
          Ab /= sqrt(Ab.rows()-1);  
        } 

        {
          // Only evaluated by 1 thread at a time
          // since it modifies global matrices
          // mutex ensures serial execution
          std::lock_guard<std::mutex> lock(pcaMutex);

          G.middleRows(start, Ab.cols()).noalias() = Ab.transpose() * Omg;
          if (i <= band / 2)
            H1.noalias() += Ab * G.middleRows(start, Ab.cols());
          else
            H2.noalias() += Ab * G.middleRows(start, Ab.cols());
          bool adjacent =
            (pi > 0 && (b + 1) == std::pow(2, pi - 1) && std::pow(2, pi) < B);
          if((!((b + 1) < band && !adjacent)) && ((i == band) || (i == band / 2) || adjacent)){
            H = H1 + H2;
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(H);
            Omg.noalias() = qr.householderQ() * Eigen::MatrixXd::Identity(n, l);
            flipOmg(Omg2, Omg);
            if (i == band) {
              H1.setZero();
              i = 0;
            } else {
              H2.setZero();
            }
          }
          i++;      
          b++;
        }
      }
    }, nchunks, min(threads, nchunks));
  }

  // get USV
  if( verbose ){
    Rcpp::Rcout << "\rFinal decompositions";
	}
  {
    Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(G);
    R.noalias() = Eigen::MatrixXd::Identity(l, m) * qr.matrixQR().triangularView<Eigen::Upper>();
    G.noalias() = qr.householderQ() * Eigen::MatrixXd::Identity(m, l);
  }
  {
    Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(G);
    Rt.noalias() = Eigen::MatrixXd::Identity(l, m) * qr.matrixQR().triangularView<Eigen::Upper>();
    G.noalias() = qr.householderQ() * Eigen::MatrixXd::Identity(m, l);
  }

  R = Rt * R;

  Eigen::MatrixXd out = R.transpose().fullPivHouseholderQr().solve(H.transpose());
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(out, Eigen::ComputeThinU | Eigen::ComputeThinV);

  if( verbose ){
    Rcpp::Rcout << "\rCompleted            " << std::endl;
  }

  Rcpp::List lst =  Rcpp::List::create(
  			Rcpp::Named("d") = Rcpp::wrap(svd.singularValues().head(k)),
        Rcpp::Named("u") = Rcpp::wrap(svd.matrixV().leftCols(k)),
        Rcpp::Named("v") = Rcpp::wrap(G * svd.matrixU().leftCols(k)));

  return lst;                            
}

