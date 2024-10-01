/***********************************************************************
 * @file		DelayedStream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		DelayedStream reads a DelayedArray into memory
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef DELAYED_STREAM_H
#define DELAYED_STREAM_H

#ifdef ARMA
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifdef EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 


#include <vector>
#include "beachmat3/beachmat.h"

using namespace Rcpp;
using namespace std;



// Handle sparse input and outputs
// Handle chunk size, workLarge.reserve(), workLarge.clear()

// returns all zeros
// res2 = GenomicDataStream:::getDA( as.matrix(M) )

/** vcfstream reads a VCF into an arma::mat in chunks, storing variants in columns.  Applies filtering for specified samples and genome region. 
 * 
*/
class DelayedStream : 
	public GenomicDataStream {
	public:

	DelayedStream(const RObject &mat) : GenomicDataStream() {

		// Get pointer to data matrix
		ptr = beachmat::read_lin_block(mat);

		// get row and col names
		mInfo = new MatrixInfo( mat );
	}


	/** destructor
	 */ 
	virtual ~DelayedStream(){
		if( mInfo != nullptr) delete mInfo;
	}


	#ifdef ARMA
	virtual bool getNextChunk( DataChunk<arma::mat, MatrixInfo> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();

		// mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem = true, strict = false)
		bool copy_aux_mem = false; // create read-only matrix without re-allocating memory
		arma::mat M(workLarge.data(), ptr->get_nrow(), ptr->get_ncol(), copy_aux_mem, true);

		chunk = DataChunk( M, *mInfo );

		return ret;
	}
	#endif

	#ifdef EIGEN
	virtual bool getNextChunk( DataChunk<Eigen::MatrixXd, 
		MatrixInfo> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(workLarge.data(), ptr->get_nrow(), ptr->get_ncol());

		chunk = DataChunk( M, *mInfo );

		return ret;
	}
	#endif

	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix, 
		MatrixInfo> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();

		Rcpp::NumericMatrix M(ptr->get_nrow(), ptr->get_ncol(), workLarge.data()); 

		chunk = DataChunk( M, *mInfo );

		return ret;
	}

	private:
	std::unique_ptr<beachmat::lin_matrix> ptr;
	vector<double> workLarge;
	bool endOfFile = false;
	MatrixInfo *mInfo = nullptr;

	bool getNextChunk_helper(){

		// set size of workspace
		vector<double> work(ptr->get_nrow());

		 // for each column j
		size_t j = 0;
	    for (j = 0; j < ptr->get_ncol(); ++j) {

	        // pointer to column i
	        auto colptr = ptr->get_col(j, work.data());
	 
	 		// insert into the larger work
	        workLarge.insert(workLarge.end(), work.begin(), work.end());
	    }
	    endOfFile = true;

	    return ! endOfFile;
	}


};











#endif