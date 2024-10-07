/***********************************************************************
 * @file		bgenstream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		reads a BGEN into matrix in chunks, storing variants in columns
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef BGEN_STREAM_H_
#define BGEN_STREAM_H_

#ifdef ARMA
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifdef EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 

#include <string>

#include <VariantInfo.h>
#include <GenomicDataStream.h>

using namespace std;

namespace GenomicDataStreamLib {

/** bgenstream reads a BGEN into an matrix in chunks, storing variants in columns.  Applies filtering for specified samples and genome region. 
 * 
*/
class bgenstream : 
	public GenomicDataStream {
	public:

	/** constructor
	*/
	bgenstream(const Param & param) : GenomicDataStream(param) {
	

	}

	/** destructor
	 */ 
	virtual ~bgenstream(){}

	#ifdef ARMA
	virtual bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk){

		return true;
	}
	#endif

	#ifdef EIGEN
	virtual bool getNextChunk( DataChunk<Eigen::MatrixXd, 
		VariantInfo> & chunk){


		return true;
	}
	#endif


	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix, 
		VariantInfo> & chunk){


		return true;
	}


	private:
	bool endOfFile = false;

	// stores genotype dosage as doubles, with the next marker inserted at the end
	// NOTE that when current size is exceeded, .insert() reallocates memory
	// this can be slow 
	// set using reserve() to set initial capacity so avaoid re-alloc
	vector<double> matDosage;	
};


}

#endif