/***********************************************************************
 * @file		GenomicDataStream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		GenomicDataStream defines an interface to read chunks of data into memory as a matrix
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef GenomicDataStream_H_
#define GenomicDataStream_H_



#ifdef ARMA
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifdef EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 

#include <boost/algorithm/string.hpp>


#include <VariantInfo.h>

using namespace std;

namespace GenomicDataStreamLib {

/** For Rcpp::NumericMatrix with MatrixInfo, set the row and col names
 */ 
static void setRowColNames( Rcpp::NumericMatrix &M, const MatrixInfo &info){
	Rcpp::rownames(M) = info.get_rownames();
	Rcpp::colnames(M) = info.get_colnames();
}



/** For all other datatypes, do nothing
 */ 
template<typename matType, typename infoType >
static void setRowColNames( matType &M, const infoType &info){
} 



/** Store the data matrix and info for a data chunk
 */ 
template<typename matType, typename infoType >
class DataChunk {
	public:
    DataChunk() : data() {}

    DataChunk( matType data) :
	    data(data) {}

    DataChunk( matType data, infoType info) :
	    data(data), info(info) { 

		// For Rcpp::NumericMatrix with MatrixInfo, 
    	// set the row and col names 
	    setRowColNames(data, info);
	 }

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






struct Param {

	Param(){}

	Param( 	const string &file,
			const string &field,
			const string &region = "",
			const string &samples = "-",
			const int &chunkSize = std::numeric_limits<int>::max(),
			const bool &missingToMean = false,
			const int &initCapacity = 200) :
		file(file), 
		field(field), 
		samples(samples), 
		chunkSize(chunkSize), 
		missingToMean(missingToMean), 
		initCapacity(initCapacity) {

		// region is string of chr:start-end delim by "\t,\n"
		// remove spaces, then split based on delim
		string region_tmp = region;
		boost::erase_all(region_tmp, " ");
    	boost::split(regions, region_tmp, boost::is_any_of("\t,\n"));

    	// remove duplicate regions
    	sort(regions.begin(), regions.end());
    	regions.erase(unique(regions.begin(), regions.end()), regions.end());

		}

	string file;
	string field;
	vector<string> regions;
	string samples;
	int chunkSize;
	bool missingToMean;
	int initCapacity;
};


/** Abstract class inheritited by vcfstream, DelayedStream
 */ 
class GenomicDataStream {
	public: 

	GenomicDataStream(){}

	/** Constructor
	 */
	GenomicDataStream( const Param & param ): param(param) {}

	/** destructor
	 */
	virtual ~GenomicDataStream() {};
	
	#ifdef ARMA
	/** Get next chunk of _features_ as arma::mat
	 * 
	 */ 
	virtual bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk){return false;
	}
	#endif

	#ifdef EIGEN
	/** Get next chunk of _features_ as Eigen::Map<Eigen::MatrixXd>
	 * 
	 */ 
	virtual bool getNextChunk( DataChunk<Eigen::Map<Eigen::MatrixXd>, VariantInfo> & chunk){return false;
	}
	#endif

	/** Get next chunk of _features_ as Rcpp::NumericMatrix
	 * 
	 */ 
	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix, VariantInfo> & chunk){return false;
	}

	/** Get next chunk of _features_ as vector<double>
	 * 
	 */ 
	virtual bool getNextChunk( DataChunk<vector<double>, VariantInfo> & chunk){return false;
	}

	protected:
	Param param;
};

}


#endif