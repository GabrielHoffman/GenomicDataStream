/***********************************************************************
 * @file		export.cpp
 * @author	   	Gabriel Hoffman
 * @email	   	gabriel.hoffman@mssm.edu
 * @brief	   	Expose GenomicDataStream library to R
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


#ifdef _OPENMP
	// [[Rcpp::plugins(openmp)]]
	#include <omp.h>
#else
	#define omp_get_num_threads() 0
	#define omp_get_thread_num() 0
#endif

#include "GenomicDataStream.h"
#include "DataTable.h"
#include "ParallelGenomicChunks.h"

#include "Rand.hpp"
#include <Rcpp/Benchmark/Timer.h>

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

	// return genotype data and variant info
	return List::create(	Named("X") = X,
							Named("info") = toDF(info) );
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
void test_DataTable(const string &file, const string &headerKey, const char delim='\t'){

	DataTable dt(file, headerKey, delim);

	dt.print(Rcpp::Rcout, "\t");
}

// [[Rcpp::export]]
Rcpp::List stream_pcaone(Rcpp::S4 gds, const string &region, int m, int k, int s = 20, int p = 7, int B = 64, int threads = 4) {
  Timer timer;
  timer.step("init params");
  // extract slots
  int n = gds.slot("nsamples"); // number of samples
  std::string file = gds.slot("file"); // fille
  std::string samples = gds.slot("samples"); // samples
  std::string field = gds.slot("field");  //  field
  int chunkSize = gds.slot("chunkSize");  //  each chunk will read max chunkSize
  double minVariance = gds.slot("minVariance");  // retain features with var > minVariance
  bool missingToMean = gds.slot("missingToMean");  
  // set region to empty here. we will used the permuted region later
  gds::Param param(file, "", samples, minVariance, chunkSize, missingToMean);
  param.setField(field);

  timer.step("init params");
  auto regions = splitRegionString( region ); // assume regions are permuted
  const size_t nchunks{regions.size()};

  MultiGenomicStreamProcessor<Eigen::MatrixXd> processor(param, regions, threads);
  // Eigen::setNbThreads(threads); 
  std::mutex pcaMutex;
  timer.step("init MultiGenomicStreamProcessor");

  const int l = k + s;
  auto randomEngine = std::default_random_engine{};
  Eigen::MatrixXd Omg = StandardNormalRandom<Eigen::MatrixXd, std::default_random_engine>(n, l, randomEngine);
  Eigen::MatrixXd Omg2 = Omg;
  Eigen::MatrixXd H1 = Eigen::MatrixXd::Zero(n, l);
  Eigen::MatrixXd H2 = Eigen::MatrixXd::Zero(n, l);
  Eigen::MatrixXd H(n, l), G(m, l), R(l, l), Rt(l, l);
  
  timer.step("init Eigen::MatrixXd for pcaone");

  size_t band = ceil((double)nchunks / B);
  
  for (int pi = 0; pi <= p; pi++) {
    if (std::pow(2, pi) >= B) {
      // reset H1, H2 to zero
      H1.setZero();
      H2.setZero();
    }
    band = std::fmin(band * 2, nchunks);

    size_t i{1},  start{0};
    processor.processChunk([&](const gds::DataChunk<Eigen::MatrixXd> &chunk, size_t b) {
      auto Ab = chunk.getData();
      standardize(Ab); // standardize
      {
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
        start += Ab.cols();
        i++;
      }
    });

    timer.step("epochs:" + to_string(pi));
  }

  // get USV
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

  timer.step("getUSV");

  Rcpp::List lst =  Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(svd.singularValues().head(k)),
                            Rcpp::Named("u") = Rcpp::wrap(svd.matrixV().leftCols(k)),
                            Rcpp::Named("v") = Rcpp::wrap(G * svd.matrixU().leftCols(k)));
  lst.attr("timing") = timer;

  return lst;                            
}
