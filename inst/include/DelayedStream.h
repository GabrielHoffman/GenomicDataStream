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
// #include "Rtatami.h"

using namespace Rcpp;
using namespace std;

namespace GenomicDataStreamLib {

// Handle sparse input and outputs
// Handle chunk size, workLarge.reserve(), workLarge.clear()

// returns all zeros
// res2 = GenomicDataStream:::getDA( as.matrix(M) )

/** Reads an Robject
 * 
*/
class DelayedStream : 
	public GenomicDataStream {
	public:

	DelayedStream( RObject mat) : GenomicDataStream() {

		Rcpp::Rcout << "read_lin_block" << std::endl;
		auto a = beachmat::read_lin_block(mat);
		Rcpp::Rcout << "end" << std::endl;


		// Get pointer to data matrix
		ptr = beachmat::read_lin_block(mat);

		NCout = ptr->get_nrow();
		NRout = ptr->get_ncol();

		// set size of workspace
		// buffer stores one row
		// buffer.reserve( ptr->get_ncol() );

		// get row and col names
		mInfo = new MatrixInfo( mat );
	}


	/** destructor
	 */ 
	~DelayedStream(){
		if( mInfo != nullptr) delete mInfo;
	}


	#ifdef ARMA
	bool getNextChunk( DataChunk<arma::mat, MatrixInfo> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();

		// mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem = true, strict = false)
		bool copy_aux_mem = false; // create read-only matrix without re-allocating memory
		arma::mat M(workLarge.data(), NRout, NCout, copy_aux_mem, true);

		chunk = DataChunk( M, *mInfo );

		return ret;
	}
	#endif

	#ifdef EIGEN
	bool getNextChunk( DataChunk<Eigen::MatrixXd, 
		MatrixInfo> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(workLarge.data(), NRout, NCout);

		chunk = DataChunk( M, *mInfo );

		return ret;
	}
	#endif

	bool getNextChunk( DataChunk<Rcpp::NumericMatrix, 
		MatrixInfo> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();


	    Rcpp::Rcout << "NumericMatrix" << std::endl;
	    Rcpp::Rcout << NRout << " " << NCout << std::endl;
	    Rcpp::Rcout << workLarge.size() << std::endl;

		Rcpp::NumericMatrix M(NRout, NCout, workLarge.data()); 


	    Rcpp::Rcout << "DataChunk" << std::endl;

		chunk = DataChunk( M, *mInfo );
	    Rcpp::Rcout << "return" << std::endl;

		return ret;
	}

	bool getNextChunk( DataChunk<vector<double>, 
		MatrixInfo> & chunk){

		// Update workLarge chunk
		bool ret = getNextChunk_helper();

		chunk = DataChunk( workLarge, *mInfo );

		return ret;
	}

	private:
	std::unique_ptr<beachmat::lin_matrix> ptr;
 
	// vector<double> buffer; 
	vector<double> workLarge;
	bool endOfFile = false;
	MatrixInfo *mInfo = nullptr;
	int NRout, NCout;

	// tatami support
	/*bool getNextChunk_helper(){

		Rcpp::Rcout << "BoundNumericPointer" << std::endl;

		Rtatami::BoundNumericPointer parsed(mat);
		Rcpp::Rcout << "success" << std::endl;
	    auto ptr = parsed->ptr;
	    NROW = ptr->nrow();
    	NCOL = ptr->ncol();

    	buffer.reserve(NCOL);

		auto wrk = ptr->dense_column();

		for (int j = 0; j < NCOL; j++) {
			auto extracted = wrk->fetch(j, buffer.data());

	 		// insert into the larger work
	        workLarge.insert(workLarge.end(), buffer.begin(), buffer.end());
		}
	    endOfFile = true;

	    return ! endOfFile;
	}*/


	// original code based on stand-alone beachmat
	bool getNextChunk_helper(){


	    vector<double> buffer(ptr->get_ncol());

		 // for each column j
		size_t j = 0;
	    for (j = 0; j < ptr->get_nrow(); j++) {


	        // pointer to column i
	        auto rp = ptr->get_row(j, buffer.data());
	 
	        for (int k =0; k<ptr->get_ncol(); k++){
	        	Rcpp::Rcout << rp[k] << ' ' << buffer[k] << std::endl;
	        }

	 		// insert into the larger work
	        workLarge.insert(workLarge.end(), buffer.begin(), buffer.end());
	    }
	    endOfFile = true;

	    Rcpp::Rcout << "end helper" << std::endl;

	    return ! endOfFile;
	}


};



}







#endif