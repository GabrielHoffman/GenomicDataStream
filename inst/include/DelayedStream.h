/***********************************************************************
 * @file		DelayedStream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		DelayedStream reads a DelayedArray into memory
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef DELAYED_STREAM_H
#define DELAYED_STREAM_H

// If this flag is not specified, run rest of code
#ifndef DISABLE_DELAYED_STREAM

#ifndef DISABLE_EIGEN
#include <Eigen/Sparse>
#endif 

#include <vector>
#include <span>

#include <Rtatami.h>

#include "GenomicDataStream_virtual.h"
#include "MatrixInfo.h"
#include "utils.h"

using namespace Rcpp;
using namespace std;

namespace gds {

/** Reads an Robject
 * 
*/
class DelayedStream : 
	public GenomicDataStream {
	public:

	DelayedStream( Rcpp::RObject robj, const vector<string> &rowNames, const int &chunkSize) 
		: GenomicDataStream(), rowNames(rowNames), chunkSize(chunkSize) {

		if( chunkSize < 1){					
			throw runtime_error("chunkSize must be positive: " + chunkSize);
		}

		parsed = new Rtatami::BoundNumericPointer(robj);

		// set current position in matrix to zero
		pos = 0;

		// set size of intermediate variables
		const auto& ptr = (*parsed)->ptr;
		NC = ptr->ncol();
		NR = ptr->nrow();

		if( ptr->nrow() != rowNames.size()){			
			throw runtime_error("DelayedStream: rowNames and nrows must be same size");
		}

		output.reserve(NC*chunkSize);
		buffer.reserve(NC);

		mInfo = new MatrixInfo();
	}


	/** destructor
	 */ 
	~DelayedStream(){
		if( mInfo != nullptr) delete mInfo;
		if( parsed != nullptr) delete parsed;
	}

	/** Get number of columns in data matrix
	 */ 
	int n_samples() override {
		return NC;
	}

	/** get FileType of param.file
	 */ 
	string getStreamType() override {
		return "DelayedStream";
	}

	bool getNextChunk( DataChunk<arma::mat> & chunk) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

	    arma::mat M(output.data(), NC, chunkSize, false, true);

		chunk = DataChunk<arma::mat>( M, mInfo );

		return ret;
	}


	bool getNextChunk( DataChunk<arma::sp_mat> & chunk) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		arma::mat M(output.data(), NC, chunkSize, false, true);

		// create sparse matrix from dense matrix
	    chunk = DataChunk<arma::sp_mat>( arma::sp_mat(M), mInfo);

		return ret;
	}

	#ifndef DISABLE_EIGEN
	bool getNextChunk( DataChunk<Eigen::MatrixXd> & chunk) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(output.data(), NC, chunkSize);

		chunk = DataChunk<Eigen::MatrixXd>( M, mInfo );

		return ret;
	}

	bool getNextChunk( DataChunk<Eigen::SparseMatrix<double> > & chunk) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(output.data(), NC, chunkSize);

		chunk = DataChunk<Eigen::SparseMatrix<double> >( M.sparseView(), mInfo );

		return ret;
	}
	#endif

	#ifndef DISABLE_RCPP
	bool getNextChunk( DataChunk<Rcpp::NumericMatrix> & chunk) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		Rcpp::NumericMatrix M(NC, chunkSize, output.data()); 
		colnames(M) = Rcpp::wrap( mInfo->getFeatureNames() );
	    // rownames(M) = Rcpp::wrap( mInfo->sampleNames );  

		chunk = DataChunk<Rcpp::NumericMatrix>( M, mInfo );

		return ret;
	}
	#endif

	bool getNextChunk( DataChunk<vector<double> > & chunk) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		chunk = DataChunk<vector<double>>( output, mInfo );

		return ret;
	}


	private:
	Rtatami::BoundNumericPointer *parsed = nullptr;	
	vector<double> buffer; 
	vector<double> output; 
	MatrixInfo *mInfo = nullptr;
	vector<string> rowNames;
	bool continueIterating = true;
	int NR, NC;
	int chunkSize;
	int pos;

	// original code based on stand-alone beachmat
	bool getNextChunk_helper(){

		// if end of file reached, return false
		if( ! continueIterating ) return continueIterating;

		// // get pointer to data
		const auto& ptr = (*parsed)->ptr;

		// // get workspace as dense row
		auto wrk = ptr->dense_row();

		// if remaning rows is less than chunkSize, 
		// 	then set chunkSize to stop at end
		chunkSize = min(chunkSize, NR - pos);

		// loop through rows
		for (int i = 0; i < chunkSize; i++) {
			// get data for row pos + i
		    auto extracted = wrk->fetch(pos + i, buffer.data());

		    // copy data into output vector in column i
		    memcpy(output.data() + NC*i, extracted, NC*sizeof(double));
		}

		// get feature names		
		mInfo->setRowNames(rowNames, pos, pos+chunkSize);

		// increment current position
		pos += chunkSize;

		// if current position is less than number of rows
		// 	return true to continue and get text chunk
		continueIterating = (pos < NR);

	    return true;
	}
};



}






#endif
#endif