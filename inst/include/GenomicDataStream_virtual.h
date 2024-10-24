/***********************************************************************
 * @file		GenomicDataStream_reader.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Define GenomicDataStream virtual class, DataChunk and Param classes
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef GENOMIC_DATA_STREAM_READER_H_
#define GENOMIC_DATA_STREAM_READER_H_

#ifndef DISABLE_EIGEN
#include <Eigen/Sparse>
#endif 

#ifndef DISABLE_RCPP
#include <RcppArmadillo.h>
#endif 

#include <boost/algorithm/string.hpp>

#include "VariantInfo.h"
#include "utils.h"

using namespace std;

namespace gds {

/** Store the data matrix and info for a data chunk
 */ 
template<typename matType, typename infoType >
class DataChunk {
	public:
    DataChunk() : data() {}

    DataChunk( matType data) :
	    data(data) {} 

    DataChunk( matType data, infoType info) :
		data(data), info(info) {}

    /** Accessor
     */ 
	matType getData() const { return data; }

	/** Accessor
     */   
	infoType getInfo() const { return info;}    

	private:
	matType data;
	infoType info;
};





/** Store parameters to pass to GenomicDataStream
 * 
 */
struct Param {

	Param(){}

	/** constructor
	 * @param file .vcf, .vcf.gz or .bcf file with tabix index
	 * @param regionString target in the format `"chr2:1-12345"`.  Multiple regions can be separated by one of `,\n\t`, for example `"chr2:1-12345, chr3:1000-8000"`. Setting region to `""` includes all variants
	 * @param samples string of comma separated sample IDs to extract: "ID1,ID2,ID3"
	 * @param chunkSize number of variants to return per chunk
	 * @param missingToMean if true, set missing values to the mean dosage value.  if false, set to NaN
	 * @param initCapacity initial capacity of temporary vector to avoid re-alloc on insert.  Size is in Mb.
	 * @param permuteFeatureOrder default is `false`. If `true` permute regions in `regionString` to avoid linkage disequilibrium betweeen nearby regions 
	 * @param rndSeed random seed for permutation
	 * 
	 * Note that field must be set separately for VCF/BCF:
	 * field: `"GT"` for genotype strings, `"DS"` for dosage, or another other field stored as an integer or float.  `"GT"` is the only string type supported
	*/
	Param( 	const string &file,
			string regionString = "",
			const string &samples = "-",
			const int &chunkSize = numeric_limits<int>::max(),
			const bool &missingToMean = false,
			const int &initCapacity = 200,
			const bool &permuteFeatureOrder = false,
			const int &rndSeed = 12345) :
		file(file), 
		samples(samples), 
		chunkSize(chunkSize), 
		missingToMean(missingToMean), 
		initCapacity(initCapacity) {

		// regionString is string of chr:start-end delim by "\t,\n"
		// remove spaces, then split based on delim
		boost::erase_all(regionString, " ");
    	boost::split(regions, regionString, boost::is_any_of("\t,\n"));

    	// remove duplicate regions, but preserve order
    	removeDuplicates( regions );

    	if( permuteFeatureOrder ){
    		// permuate order of region to avoid correlated features
    		// in streaming SVD
    		std::shuffle( regions.begin(), regions.end(), std::mt19937(rndSeed));
    	}
	}

	void setField( const string &field_){
		field = field_;
	}

	string file;
	string field = "";
	vector<string> regions;
	string samples;
	int chunkSize;
	bool missingToMean;
	int initCapacity;
};


/** Virtual class inheritited by vcfstream, bgenstream, DelayedStream
 */ 
class GenomicDataStream {
	public: 

	GenomicDataStream() {}

	/** Constructor
	 */
	GenomicDataStream( const Param & param ) : param(param) {}

	/** destructor
	 */
	virtual ~GenomicDataStream() {};
	
	/** Get number of columns in data matrix
	 */ 
	virtual int n_samples() = 0;

	/** Get next chunk of _features_ as arma::mat
	 */ 
	virtual bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk) = 0;

	/** Get next chunk of _features_ as arma::sp_mat
	 */ 
	virtual bool getNextChunk( DataChunk<arma::sp_mat, VariantInfo> & chunk) = 0;

	#ifndef DISABLE_EIGEN
	/** Get next chunk of _features_ as Eigen::MatrixXd
	 */ 
	virtual bool getNextChunk( DataChunk<Eigen::MatrixXd, VariantInfo> & chunk) = 0;

	/** Get next chunk of _features_ as SparseMatrix<double>
	 */ 
	virtual bool getNextChunk( DataChunk<Eigen::SparseMatrix<double>, VariantInfo> & chunk) = 0;
	#endif

	#ifndef DISABLE_RCPP
	/** Get next chunk of _features_ as Rcpp::NumericMatrix
	 */ 
	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix, VariantInfo> & chunk)  = 0;
	#endif

	/** Get next chunk of _features_ as vector<double>
	 * 
	 */ 
	virtual bool getNextChunk( DataChunk<vector<double>, VariantInfo> & chunk) = 0;

	protected:
	Param param;
};



}


#endif