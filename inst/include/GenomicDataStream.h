/***********************************************************************
 * @file		GenomicDataStream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		GenomicDataStream defines an interface to read chunks of data into memory as a matrix
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

/*! \mainpage A scalable interface between data and analysis underneath R
 *
 \image html "GenomicDataStream.png"
<div style="text-align: justify">
Reading genomic data files ([VCF](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/), [BCF](https://samtools.github.io/bcftools/howtos/index.html), [BGEN](https://www.chg.ox.ac.uk/~gav/bgen_format/index.html), [H5AD](https://anndata.readthedocs.io/en/latest/index.html), [DelayedArray](https://bioconductor.org/packages/DelayedArray)) into R/Rcpp in chunks for analysis with [Armadillo](https://doi.org/10.21105/joss.00026) / [Eigen](eigen.tuxfamily.org) / [Rcpp](https://www.rcpp.org) libraries.  Mondern datasets are often too big to fit into memory, and many analyses operate a small chunk features at a time.  Yet in practice, many implementations require the whole dataset stored in memory.  Others pair an analysis with a specific data format (i.e. regresson analysis paired with genotype data from a VCF) in way that the two components can't be separated for use in other applications. 
 </div>
 The `GenomicDataStream` C++ inferface separates 

 -# data source 
 -# streaming chunks of features into a data matrix
 -# downstream analysis  
 * 
 * 
 * ## Dependencies
| Package | Ref | Role |
| - | --- | --------- |
[vcfppR](https://cran.r-project.org/package=vcfppR) | [Bioinformatics](https://doi.org/10.1093/bioinformatics/btae049)  | C++ API for htslib  |
[htslib](https://github.com/samtools/htslib) | [GigaScience](https://doi.org/10.1093/gigascience/giab007)  | C API for VCF/BCF files |
[beatchmat](https://bioconductor.org/packages/beachmat/) | [PLoS Comp Biol](https://doi.org/10.1371/journal.pcbi.1006135)  | C++ API for access data owned by R |
[Rcpp](https://cran.r-project.org/package=Rcpp)| [J Stat Software](https://doi.org/10.18637/jss.v040.i08) |  API for R/C++ integration
[RcppEigen](https://cran.r-project.org/package=RcppEigen) | [J Stat Software](https://doi.org/10.18637/jss.v052.i05) | API for Rcpp access to Eigen matrix library
[RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo)| [J Stat Software](https://doi.org/10.18637/jss.v040.i08) | API for Rcpp access to Armadillo matrix library
[Eigen](eigen.tuxfamily.org) | |C++ library for linear algebra with advanced features
[Armadillo](https://arma.sourceforge.net) | [J Open Src Soft](https://doi.org/10.21105/joss.00026) | User-friendly C++ library for linear algebra


 * 
 * 
 */

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

#include "VariantInfo.h"
#include "utils.h"

using namespace std;

namespace GenomicDataStreamLib {


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





/** Store parameters to pass to GenomicDataStream
 * 
 */
struct Param {

	Param(){}

	/** constructor
	 * @param file .vcf, .vcf.gz or .bcf file with tabix index
	 * @param field `"GT"` for genotype strings, `"DS"` for dosage, or another other field stored as an integer or float.  `"GT"` is the only string type supported
	 * @param regionString target in the format `"chr2:1-12345"`.  Multiple regions can be separated by one of `,\n\t`, for example `"chr2:1-12345, chr3:1000-8000"`. Setting region to `""` includes all variants
	 * @param samples string of comma separated sample IDs to extract: "ID1,ID2,ID3"
	 * @param chunkSize number of variants to return per chunk
	 * @param missingToMean if true, set missing values to the mean dosage value.  if false, set to NaN
	 * @param initCapacity initial capacity of temporary vector to avoid re-alloc on insert.  Size is in Mb.
	 * @param permuteFeatureOrder default is `false`. If `true` permute regions in `regionString` to avoid linkage disequilibrium betweeen nearby regions 
	 * @param rndSeed random seed for permutation
	*/
	Param( 	const string &file,
			const string &field,
			string regionString = "",
			const string &samples = "-",
			const int &chunkSize = std::numeric_limits<int>::max(),
			const bool &missingToMean = false,
			const int &initCapacity = 200,
			const bool &permuteFeatureOrder = false,
			const int &rndSeed = 12345) :
		file(file), 
		field(field), 
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

	string file;
	string field;
	vector<string> regions;
	string samples;
	int chunkSize;
	bool missingToMean;
	int initCapacity;
};


/** Abstract class inheritited by vcfstream, bgenstream, DelayedStream
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
	
	/** Get number of columns in data matrix
	 */ 
	virtual int n_cols(){};

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