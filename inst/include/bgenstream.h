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

		View::UniquePtr bgenView;
		bgenView = View::create( param.file );

		// define range query
		// view->set_query( query ) ;


		std::size_t const number_of_variants = bgenView->number_of_variants();
		// std::size_t const number_of_samples = bgenView->number_of_samples();
		std::size_t number_of_samples;

		std::string SNPID, rsid, chromosome;
		genfile::bgen::uint32_t position;
		std::vector< std::string > alleles;

		int max_entries_per_sample = 1e7;
		Dimension data_dimension = Dimension( number_of_variants, number_of_samples, max_entries_per_sample ) ;
		Dimension ploidy_dimension = Dimension( number_of_variants, number_of_samples ) ;

		NumericVector data = NumericVector( data_dimension, NA_REAL ) ;
		IntegerVector ploidy = IntegerVector( ploidy_dimension, NA_INTEGER ) ;
		LogicalVector phased = LogicalVector( number_of_variants, NA_LOGICAL ) ;

		std::vector< std::string > sampleNames ;
		std::map< std::size_t, std::size_t > requestedSamplesByIndexInDataIndex ;
	
		// if( requestedSamples ) {
		// 	get_requested_samples( *view, *requestedSamples, &number_of_samples, &sampleNames, &requestedSamplesByIndexInDataIndex ) ;
		// } else {
			get_all_samples( *bgenView, &number_of_samples, &sampleNames, &requestedSamplesByIndexInDataIndex ) ;
		// }

		// Iterate through variants
		for( std::size_t variant_i = 0; variant_i < number_of_variants; ++variant_i ) {
			bgenView->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ;

			Rcpp::Rcout << SNPID << std::endl;
			// chromosomes[variant_i] = chromosome ;
			// positions[variant_i] = position ;
			// rsids[variant_i] = rsid ;
			// number_of_alleles[variant_i] = alleles.size() ;
			// allele0s[variant_i] = alleles[0] ;
			// allele1s[variant_i] = alleles[1] ;

			// DataSetter setter(
			// 	&ploidy, ploidy_dimension,
			// 	&data, data_dimension,
			// 	&phased,
			// 	variant_i,
			// 	requestedSamplesByIndexInDataIndex
			// ) ;
			
			// bgenView->read_genotype_data_block( setter ) ; // will be fixed later
		}

	}

	/** destructor
	 */ 
	virtual ~bgenstream(){}

	#ifdef ARMA
	virtual bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk){

		return true;
	}
	#endif

	#ifdef EIGEN
	virtual bool getNextChunk( DataChunk<Eigen::MatrixXd, 
		VariantInfo> & chunk){


		return true;
	}
	#endif


	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix, 
		VariantInfo> & chunk){


		return true;
	}


	private:
	bool endOfFile = false;

	// stores genotype dosage as doubles, with the next marker inserted at the end
	// NOTE that when current size is exceeded, .insert() reallocates memory
	// this can be slow 
	// set using reserve() to set initial capacity so avaoid re-alloc
	vector<double> matDosage;	
};


}

#endif