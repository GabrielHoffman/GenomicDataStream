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



/** vcfstream reads a VCF into an arma::mat in chunks, storing variants in columns.  Applies filtering for specified samples and genome region. 
 * 
*/
template<typename T>
class DelayedStream : 
	public GenomicDataStream {
	public:

	DelayedStream(const RObject &mat) : GenomicDataStream() {

		T ptr = beachmat::read_lin_block(mat);


		// vector<double> workspace(ptrTmp->get_nrow());

	}


	/** destructor
	 */ 
	virtual ~DelayedStream(){
	}


	#ifdef ARMA
	virtual bool getNextChunk( DataChunk<arma::mat, int> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();

		// mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem = true, strict = false)
		bool copy_aux_mem = false; // create read-only matrix without re-allocating memory
		arma::mat M(workLarge.data(), ptr->get_nrow(), ptr->get_ncol(), copy_aux_mem, true);

		chunk = DataChunk<arma::mat, int>( M );

		return ret;
	}
	#endif

	#ifdef EIGEN
	virtual bool getNextChunk( DataChunk<Eigen::Map<Eigen::MatrixXd>, 
		int> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();

		Eigen::Map<Eigen::MatrixXd> M(workLarge.data(), ptr->get_nrow(), ptr->get_ncol());

		chunk = DataChunk<Eigen::Map<Eigen::MatrixXd>, int>( M );

		return ret;
	}
	#endif

	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix, 
		int> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();

		Rcpp::NumericMatrix M(ptr->get_nrow(), ptr->get_ncol(), workLarge.data()); 
		chunk = DataChunk<Rcpp::NumericMatrix, int>( M );

		return ret;
	}

	private:
	T ptr;
	vector<double> workspace, workLarge;

	bool getNextChunk_helper(){
		 // for each column i
	    for (size_t i = 0; i < ptr->get_ncol(); ++i) {

	        // pointer to column i
	        auto colptr = ptr->get_col(i, workspace.data());
	 
	 		// insert into the larger workspace
	        workLarge.insert(workspace);
	    }

	    return true;
	}


};













#endif