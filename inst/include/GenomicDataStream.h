/***********************************************************************
 * @file		GenomicDataStream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		GenomicDataStream defines an interface to read chunks of data into memory
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef GenomicDataStream_H_
#define GenomicDataStream_H_


#include <RcppArmadillo.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]


#include <VariantInfo.h>

using namespace std;


template<typename matType, typename infoType >
class DataChunk {
	public:
    DataChunk() {}

    DataChunk( matType data, infoType info) :
	    data(data), info(info) {}

	matType data;
	infoType info;
};




struct Param {

	Param( 	const string &file,
			const string &field,
			const string &region = "",
			const string &samples = "-",
			const int &chunkSize = std::numeric_limits<int>::max(),
			const bool &missingToMean = false,
			const int &initCapacity = 200) :
		file(file), 
		field(field), 
		region(region), 
		samples(samples), 
		chunkSize(chunkSize), 
		missingToMean(missingToMean), 
		initCapacity(initCapacity) {}

	string file;
	string field;
	string region;
	string samples;
	int chunkSize;
	bool missingToMean;
	int initCapacity;
};


class GenomicDataStream {
	public: 

	/** Constructor
	 */
	GenomicDataStream( const Param & param ): param(param) {}

	/** destructor
	 */
	virtual ~GenomicDataStream(){}
	
	/** Get next chunk of _features_ as arma::mat
	 * 
	 */ 
	virtual bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk){return false;
	}

	/** Get next chunk of _features_ as Eigen::Map<Eigen::MatrixXd>
	 * 
	 */ 
	virtual bool getNextChunk( DataChunk<Eigen::Map<Eigen::MatrixXd>, VariantInfo> & chunk){return false;
	}

	protected:
	Param param;
};




#endif