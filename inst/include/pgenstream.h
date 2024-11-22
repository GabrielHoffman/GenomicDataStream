/***********************************************************************
 * @file		pgenstream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		reads a plink2/PGEN into matrix in chunks, storing variants in columns
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef PGEN_STREAM_H_
#define PGEN_STREAM_H_

#ifndef DISABLE_EIGEN
#include <Eigen/Sparse>
#endif 

#include <string>

#include "VariantInfo.h"
#include "GenomicDataStream_virtual.h"
#include "GenomicRanges.h"

// #include "pgen/RPgenReader.h"
#include "pgen/pvar.h"


// RPgenReader *pg = nullptr;
// RPvar *pvar = nullptr;

using namespace std;
using namespace arma;

namespace gds {

/** pgenstream reads a PGEN into an matrix in chunks, storing variants in columns.  Applies filtering for specified samples and genome region. 
 * 
*/
class pgenstream : 
	public GenomicDataStream {
	public:

	pgenstream() {}

	/** constructor
	*/
	pgenstream(const Param & param) : GenomicDataStream(param) {
		
		// Read index file (pvar/bim)
		// pvar = new RPvar();
        // pvar->Load(param.fileIdx);


		// pg = new RPgenReader();

		// Nullable<List> pv = R_NilValue;
		// Nullable<int> raw_sample_ct = R_NilValue;
        // Nullable<IntegerVector> sample_subset = R_NilValue;

		// // Read data file (pgen/bed)
		// pg = new RPgenReader();
  		// pg->Load(param.file, pv, raw_sample_ct, sample_subset);


  		// // read data into matrix
  		// NumericMatrix buf;
  		// IntegerVector varIdx;
  		// pg->ReadList( buf, varIdx, param.missingToMean);





  		// // Change
  		// // RPgenReader::ReadList() to take vector<> instead of R types
  		// // Does this handle imputed data? or just int?
  		// // convert GenomicRanges to idx
  		// pg->ReadList( dsg, idx, param.missingToMean);
	


	}

	/** destructor
	 */ 
	~pgenstream(){
		if( vInfo != nullptr) delete vInfo;
		// if( pg != nullptr){
		// 	pg->Close();
		// 	delete pg;
		// }
		// if( pvar != nullptr){
		// 	pvar->Close();
		// 	delete pvar;
		// }
	}

	/** Get number of columns in data matrix
	 */ 
	int n_samples() override {
		return number_of_samples;
	}

	/** get FileType of param.file
	 */ 
	string getStreamType() override {
		return toString( param.fileType);
	}

	bool getNextChunk( DataChunk<arma::mat> & chunk) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		arma::mat M(matDosage.data(), number_of_samples, vInfo->size(), false, true);

	    chunk = DataChunk<arma::mat>( M, vInfo );

		return ret;
	}

	bool getNextChunk( DataChunk<arma::sp_mat> & chunk) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		arma::mat M(matDosage.data(), number_of_samples, vInfo->size(), false, true);

	    chunk = DataChunk<arma::sp_mat>( arma::sp_mat(M), vInfo );

		return ret;
	}

	#ifndef DISABLE_EIGEN
	bool getNextChunk( DataChunk<Eigen::MatrixXd> & chunk) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(matDosage.data(), number_of_samples, vInfo->size());

		chunk = DataChunk<Eigen::MatrixXd>( M, vInfo );

		return ret;
	}


	bool getNextChunk( DataChunk<Eigen::SparseMatrix<double> > & chunk) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(matDosage.data(), number_of_samples, vInfo->size());

		chunk = DataChunk<Eigen::SparseMatrix<double>>( M.sparseView(), vInfo );

		return ret;
	}
	#endif

	#ifndef DISABLE_RCPP
	bool getNextChunk( DataChunk<Rcpp::NumericMatrix> & chunk) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Rcpp::NumericMatrix M(number_of_samples, vInfo->size(), matDosage.data()); 
		colnames(M) = Rcpp::wrap( vInfo->getFeatureNames() );
	    rownames(M) = Rcpp::wrap( vInfo->sampleNames );  

		chunk = DataChunk<Rcpp::NumericMatrix>( M, vInfo );

		return ret;
	}
	#endif

	bool getNextChunk( DataChunk<vector<double> > & chunk) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		chunk = DataChunk<vector<double> >( matDosage, vInfo );

		return ret;
	}

	private:
	size_t number_of_samples = 0;
	vector<double> matDosage;
	VariantInfo *vInfo = nullptr;
	// RPgenReader *pg = nullptr;
	RPvar *pvar = nullptr;

	bool getNextChunk_helper(){	


  		// pg->ReadList( matDosage, idx, param.missingToMean);
	

		return true;
	}
};

}

#endif