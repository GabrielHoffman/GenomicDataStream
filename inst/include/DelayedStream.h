/***********************************************************************
 * @file		DelayedStream.h
 * @author	Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		DelayedStream reads a DelayedArray into memory in chunks
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
#include "beachmat3/beachmat.h"
#include "tatami/tatami.hpp"

#include "GenomicDataStream_virtual.h"
#include "MatrixInfo.h"
#include "utils.h"
#include "GenomicRanges.h"

using namespace Rcpp;
using namespace std;

namespace gds {

/** Reads an Robject
 * 
*/
class DelayedStream : 
	public GenomicDataStream {
	public:

	DelayedStream( 
		const shared_ptr<tatami::NumericMatrix> &ptr, 
		const vector<string> &rowNames, 
		const int &chunkSize) 
		: GenomicDataStream(), 
		ptr(ptr), 
		rowNames(rowNames), chunkSize(chunkSize) 
		{

		if( chunkSize < 1){					
			throw runtime_error("chunkSize must be positive: " + to_string(chunkSize));
		}
		// set current position in matrix to zero
		pos = 0;

		// set size of intermediate variables
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
	}

	/** setter
	 */
	void setRegions(const vector<string> &regions) override {
		throw runtime_error("DelayedStream: setRegions() not implemented");
	}

	/** Get number of columns in data matrix
	 */ 
	int n_samples() override {
		return NC;
	}

	/** Get number of rows in data matrix
	 */ 
	int n_rows() {
		return NR;
	}

	/** Get vector of sample names in order that the genotypes are extracted
	 */ 
	vector<string> getSampleNames() override {
		throw runtime_error("DelayedStream: getSampleNames() not implemented");
		return vector<string>(1);
	}

	/** get FileType of param.file
	 */ 
	string getStreamType() override {
		return "DelayedStream";
	}

	GenomicRanges getChromRanges() override {
		return GenomicRanges();
	}

	bool getNextChunk( DataChunk<arma::mat> & chunk, const bool &useFilter = true) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		arma::mat M(output.data(), NC, chunkSize, false, true);

		chunk = DataChunk<arma::mat>( M, mInfo );

		return ret;
	}

	bool getNextChunk( DataChunk<arma::sp_mat> & chunk, const bool &useFilter = true) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		arma::mat M(output.data(), NC, chunkSize, false, true);

		// create sparse matrix from dense matrix
		chunk = DataChunk<arma::sp_mat>( arma::sp_mat(M), mInfo);

		return ret;
	}

	#ifndef DISABLE_EIGEN
	bool getNextChunk( DataChunk<Eigen::MatrixXd> & chunk, const bool &useFilter = true) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(output.data(), NC, chunkSize);

		chunk = DataChunk<Eigen::MatrixXd>( M, mInfo );

		return ret;
	}

	bool getNextChunk( DataChunk<Eigen::SparseMatrix<double> > & chunk, const bool &useFilter = true) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(output.data(), NC, chunkSize);

		chunk = DataChunk<Eigen::SparseMatrix<double> >( M.sparseView(), mInfo );

		return ret;
	}
	#endif

	#ifndef DISABLE_RCPP
	bool getNextChunk( DataChunk<Rcpp::NumericMatrix> & chunk, const bool &useFilter = true) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		Rcpp::NumericMatrix M(NC, chunkSize, output.data()); 
		colnames(M) = Rcpp::wrap( mInfo->getFeatureNames() );
    // rownames(M) = Rcpp::wrap( mInfo->sampleNames );  

		chunk = DataChunk<Rcpp::NumericMatrix>( M, mInfo );

		return ret;
	}
	#endif

	bool getNextChunk( DataChunk<vector<double> > & chunk, const bool &useFilter = true) override {

		// Update vector<double> output
		bool ret = getNextChunk_helper();

		chunk = DataChunk<vector<double>>( output, mInfo );

		return ret;
	}

	/**
	 * Get chunks based on start and end 
	 */ 
	void getNextChunk( DataChunk<Eigen::MatrixXd> & chunk, int start, int len){

		int length = min(len, NR - start);

		// get consecutive rows
		// false = row extractor
    auto wrk = tatami::consecutive_extractor<false>(ptr.get(), true, pos, chunkSize);

		// loop through rows
		for (int i = 0; i < length; i++) {
			// using consecutive_extractor, don't need index
	    auto extracted = wrk->fetch(buffer.data());

	    // copy data into output vector in column i
	    memcpy(output.data() + NC*i, extracted, NC*sizeof(double));
		}

		// get feature names		
		mInfo->setRowNames(rowNames, start, start + length);

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(output.data(), NC, length);

		chunk = DataChunk<Eigen::MatrixXd>( M, mInfo );
	}

	private:
	const shared_ptr<tatami::NumericMatrix> &ptr;
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

		// get workspace as dense row
		// auto wrk = ptr->dense_row();

		// get consecutive rows
		// false = row extractor
    auto wrk = tatami::consecutive_extractor<false>(ptr.get(), true, pos, chunkSize);

		// if remaning rows is less than chunkSize, 
		// 	then set chunkSize to stop at end
		chunkSize = min(chunkSize, NR - pos);

		// loop through rows
		for (int i = 0; i < chunkSize; i++) {
			// get data for row pos + i
	    // auto extracted = wrk->fetch(pos + i, buffer.data());

			// using consecutive_extractor, don't need index
	    auto extracted = wrk->fetch(buffer.data());

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
