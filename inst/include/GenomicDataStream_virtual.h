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

#include <string>
#include <filesystem>

#include "VariantInfo.h"
#include "utils.h"

using namespace std;

namespace gds {

/** Store the data matrix and info for a data chunk
 */ 
template<typename matType>
class DataChunk {
	public:
    DataChunk() : data() {}

    DataChunk( matType data) :
	    data(data) {} 

    DataChunk( matType data, DataInfo *info) :
		data(data), info(info) {}

    /** Accessor
     */ 
	matType getData() const { return data; }

	/** Accessor
     */   
	template<typename infoType>
	infoType * getInfo() { 
		return static_cast<infoType *>( info );
	}    

	// private:
	matType data;
	DataInfo *info;
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
	 * @param minVariance features with variance >= minVariance are retained
	 * @param permuteFeatureOrder default is `false`. If `true` permute regions in `regionString` to avoid linkage disequilibrium betweeen nearby regions 
	 * @param rndSeed random seed for permutation
	 * 
	 * Note that field must be set separately for VCF/BCF:
	 * field: `"GT"` for genotype strings, `"DS"` for dosage, or another other field stored as an integer or float.  `"GT"` is the only string type supported
	*/
	Param( 	const string &file,
			string regionString = "",
			const string &samples = "-",
			const double minVariance = 0,
			const int &chunkSize = numeric_limits<int>::max(),
			const bool &missingToMean = false,
			const int &initCapacity = 200,
			const bool &permuteFeatureOrder = false,
			const int &rndSeed = 12345) :
		file( std::filesystem::absolute(file) ), 
		samples(samples), 
		chunkSize(chunkSize), 
		missingToMean(missingToMean), 
		initCapacity(initCapacity),
		minVariance(minVariance),
		fileType(getFileType(file)) {

		if( chunkSize < 1){					
			throw runtime_error("chunkSize must be positive: " + to_string(chunkSize));
		}

		if( minVariance < 0){
			throw runtime_error("minVariance must be positive");
		}

		// parse regions
		regions = splitRegionString( regionString );

    	if( permuteFeatureOrder ){
    		// permuate order of region to avoid correlated features
    		// in streaming SVD
    		std::shuffle( regions.begin(), regions.end(), std::mt19937(rndSeed));
    	}
	}

	void setField( const string &field_) {
		field = field_;
	}

	/** Custom path to PSAM/FAM file
	 */
	void setSamplesFile( const string &file) {
		fileSamples = file;
	}

	string file, fileSamples = "";
	string field = "";
	vector<string> regions;
	string samples;
	int chunkSize;
	bool missingToMean;
	int initCapacity;
	double minVariance;
	FileType fileType;
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

	/** setter
	 */
	virtual void setRegions(const vector<string> &regions) = 0; 
	
	/** Get number of columns in data matrix
	 */ 
	virtual int n_samples() = 0;

	/** get FileType of param.file
	 */ 
	virtual string getStreamType() = 0;

	/** get minVariance stored in Param
	 */ 
	double getMinVariance(){
		return param.minVariance;
	}

	/** Get next chunk of _features_ as arma::mat
	 */ 
	virtual bool getNextChunk( DataChunk<arma::mat> & chunk) = 0;

	/** Get next chunk of _features_ as arma::sp_mat
	 */ 
	virtual bool getNextChunk( DataChunk<arma::sp_mat> & chunk) = 0;

	#ifndef DISABLE_EIGEN
	/** Get next chunk of _features_ as Eigen::MatrixXd
	 */ 
	virtual bool getNextChunk( DataChunk<Eigen::MatrixXd> & chunk) = 0;

	/** Get next chunk of _features_ as SparseMatrix<double>
	 */ 
	virtual bool getNextChunk( DataChunk<Eigen::SparseMatrix<double> > & chunk) = 0;
	#endif

	#ifndef DISABLE_RCPP
	/** Get next chunk of _features_ as Rcpp::NumericMatrix
	 */ 
	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix> & chunk)  = 0;
	#endif

	/** Get next chunk of _features_ as vector<double>
	 * 
	 */ 
	virtual bool getNextChunk( DataChunk<vector<double> > & chunk) = 0;

	protected:
	Param param;
};

/** Given matDosage.data(), number_of_samples, vInfo->size()
 * set matDosage and vInfo so the first K entries are the valid
 * features
*/
static void applyVarianceFilter(vector<double> &matDosage, VariantInfo *vInfo, const int &number_of_samples, const double &minVariance = 0){

	if( minVariance > 0 ){
		// create temp arma matrix
		arma::mat M(matDosage.data(), number_of_samples, vInfo->size(), false, true);

		// variance by columns
		auto colvars = var(M);

		// indeces passing minimum variance cutoff
		arma::uvec idx = find(colvars >= minVariance);

		// convert uvec to vector<unsigned int>
		vector<unsigned int> idx2(idx.n_elem);
	    copy(idx.begin(), idx.end(), idx2.begin());

		// keep only variants specified in idx2
		vInfo->retainVariants( idx2 );

		// set subset of columns
		arma::mat M_subset = M.cols(idx);

		// copy M_subset into matDosage
		memcpy(matDosage.data(), M_subset.memptr(), M_subset.n_elem*sizeof(double));
	}
}




} // end namespace
#endif
