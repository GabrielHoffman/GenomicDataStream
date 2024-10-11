/***********************************************************************
 * @file		bgenstream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		reads a BGEN into matrix in chunks, storing variants in columns
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef BGEN_STREAM_H_
#define BGEN_STREAM_H_

#ifdef ARMA
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifdef EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 

#include <string>

#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"

#include <boost/algorithm/string.hpp>

#include <VariantInfo.h>
#include <GenomicDataStream.h>
#include "load.h"

using namespace std;
using namespace genfile::bgen;
using namespace Rcpp;

namespace GenomicDataStreamLib {

/** bgenstream reads a BGEN into an matrix in chunks, storing variants in columns.  Applies filtering for specified samples and genome region. 
 * 
*/
class bgenstream : 
	public GenomicDataStream {
	public:

	/** constructor
	*/
	bgenstream(const Param & param) : GenomicDataStream(param) {

		Rcpp::DataFrame const ranges;
		vector<string> const requested_rsids;
		size_t max_entries_per_sample;	

		view = construct_view( param.file, param.file + ".bgi", ranges, requested_rsids ) ;

		number_of_variants = view->number_of_variants() ;

		// if get all samples
		if( param.samples.compare("-") == 0 ){
			get_all_samples( *view, &number_of_samples, &sampleNames, &requestedSamplesByIndexInDataIndex ) ;
		}else{

			// get subset of samples
			vector<string> requestedSamples;

			// boost::erase_all(param.samples, " ");
    		boost::split(requestedSamples, param.samples, boost::is_any_of("\t,\n"));

			get_requested_samples( *view, requestedSamples, &number_of_samples, &sampleNames, &requestedSamplesByIndexInDataIndex ) ;
		}


		vInfo = new VariantInfo( sampleNames );

		matDosage.reserve( number_of_variants*number_of_samples );

	}

	/** destructor
	 */ 
	~bgenstream(){
		if( vInfo != nullptr) delete vInfo;
		if( data != nullptr) delete data;
	}

	#ifdef ARMA
	virtual bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk){


		Rcpp::Rcout << "getNextChunk_helper" << endl;
		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Rcpp::Rcout << "matDosage" << endl;
		for(size_t i=0; i<data->size(); i++){
			if(i < 10) Rcpp::Rcout << data->at(i) << ' ';
			matDosage[i] = data->at(i);
		}

		Rcpp::Rcout << "\nInit matrix" << endl;
		// mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem = true, strict = false)
		bool copy_aux_mem = false; // create read-only matrix without re-allocating memory
		arma::mat M(matDosage.data(), number_of_samples, number_of_variants, copy_aux_mem, true);

		Rcpp::Rcout << "DataChunk" << endl;
	    chunk = DataChunk<arma::mat, VariantInfo>( M, *vInfo );

		return true;
	}
	#endif

	#ifdef EIGEN
	virtual bool getNextChunk( DataChunk<Eigen::MatrixXd, VariantInfo> & chunk){


		return true;
	}
	#endif


	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix, VariantInfo> & chunk){


		return true;
	}


	virtual bool getNextChunk( DataChunk<vector<double>, VariantInfo> & chunk){


		return true;
	}




	private:
	View::UniquePtr view; 
	std::size_t number_of_variants;
	std::size_t number_of_samples = 0;
	std::vector< std::string > sampleNames;
	std::map< std::size_t, std::size_t > requestedSamplesByIndexInDataIndex;
	VariantInfo *vInfo = nullptr;
	vector<double> *data = nullptr;
		

	// stores genotype dosage as doubles, with the next marker inserted at the end
	// NOTE that when current size is exceeded, .insert() reallocates memory
	// this can be slow 
	// set using reserve() to set initial capacity so avaoid re-alloc
	vector<double> matDosage;


	bool getNextChunk_helper(){	

		// Declare storage for all the things we need
		StringVector chromosomes( number_of_variants ) ;
		IntegerVector positions( number_of_variants ) ;
		StringVector rsids( number_of_variants ) ;
		IntegerVector number_of_alleles( number_of_variants ) ;
		StringVector allele0s( number_of_variants ) ;
		StringVector allele1s( number_of_variants ) ;

		int max_entries_per_sample = 3;
		Dimension data_dimension = Dimension( number_of_variants, number_of_samples, max_entries_per_sample ) ;
		Dimension ploidy_dimension = Dimension( number_of_variants, number_of_samples ) ;

		// NumericVector data = NumericVector( data_dimension, NA_REAL ) ;
		data = new vector<double> ( number_of_variants * number_of_samples * max_entries_per_sample);
		IntegerVector ploidy = IntegerVector( ploidy_dimension, NA_INTEGER ) ;
		LogicalVector phased = LogicalVector( number_of_variants, NA_LOGICAL ) ;

		string SNPID, rsid, chromosome;
		genfile::bgen::uint32_t position;
		vector<string> alleles;

		// Iterate through variants
		for( size_t j = 0; j < number_of_variants; j++ ) {

			Rcpp::Rcout << j << endl;

			view->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ;
			chromosomes[j] = chromosome ;
			positions[j] = position ;
			rsids[j] = rsid ;
			number_of_alleles[j] = alleles.size() ;
			allele0s[j] = alleles[0] ;
			allele1s[j] = alleles[1] ;

			DataSetter setter(
				&ploidy, ploidy_dimension,
				data, data_dimension,
				&phased,
				j,
				requestedSamplesByIndexInDataIndex
			) ;
			
			view->read_genotype_data_block( setter ) ;		
		}



		DataFrame variants = DataFrame::create(
			Named("chromosome") = chromosomes,
			Named("position") = positions,
			Named("rsid") = rsids,
			Named("number_of_alleles") = number_of_alleles,
			Named("allele0") = allele0s,
			Named("allele1") = allele1s
		) ;
		variants.attr( "row.names" ) = rsids ;

		StringVector genotypeNames(max_entries_per_sample) ;
		for( std::size_t i = 0; i < max_entries_per_sample; ++i ) {
			genotypeNames[i] = "g=" + atoi(i) ;
		}

		List dataNames = List(3) ;
		dataNames[0] = rsids ;
		dataNames[1] = sampleNames ;
		dataNames[2] = genotypeNames ;

		List ploidyNames = List(2) ;
		ploidyNames[0] = rsids ;
		ploidyNames[1] = sampleNames ;

		// data.attr( "dimnames" ) = dataNames ;
		ploidy.attr( "dimnames" ) = ploidyNames ;

		List result ;
		result[ "variants" ] = variants ;
		result[ "samples" ] = sampleNames ;
		result[ "ploidy" ] = ploidy ;
		result[ "phased" ] = phased ;
		// result[ "data" ] = data ;

		// matDosage = vector<double>(data);

		// return( result ) ;
		return false;
	}
};


}

#endif