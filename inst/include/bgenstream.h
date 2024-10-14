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


/** TODO
 * chunks
 * Remove Rcpp dependency
 * Currently uses Armadillo cube
 * How to convert to Eigen, NumericMatrix, vector<double>
 * 
 */

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

		// Initialize view
		view = construct_view( param.file, param.file + ".bgi", ranges, requested_rsids ) ;

		// Filter variants
		IndexQuery::UniquePtr query = IndexQuery::create( param.file + ".bgi" );
		vector<string> rsids = {"RSID_101", "RSID_2", "RSID_102"};
		query->include_rsids( rsids ) ;
		query->initialise() ;

		view->set_query( query ) ;

		number_of_variants = view->number_of_variants() ;

		// Filter samples
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

		// store probabilities
		int n = 1e6 * param.initCapacity / (double) (sizeof(double) * number_of_samples * max_entries_per_sample);
		probs.reserve( n );

		// store dosage
		n = 1e6 * param.initCapacity / (double) (sizeof(double) * number_of_samples);
		matDosage.reserve(n);
	}

	/** destructor
	 */ 
	~bgenstream(){
		if( vInfo != nullptr) delete vInfo;
	}

	#ifdef ARMA
	virtual bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk){

		Rcpp::Rcout << "getNextChunk_helper()..." << endl;
		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		arma::mat M(matDosage.data(), number_of_samples, vInfo->size(), false, true);

	    chunk = DataChunk<arma::mat, VariantInfo>( M, *vInfo );


		return ret;
	}
	#endif

	#ifdef EIGEN
	virtual bool getNextChunk( DataChunk<Eigen::MatrixXd, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(matDosage.data(), number_of_samples, vInfo->size());

		chunk = DataChunk<Eigen::MatrixXd, VariantInfo>( M, *vInfo );

		return ret;
	}
	#endif


	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Rcpp::NumericMatrix M(number_of_samples, vInfo->size(), matDosage.data()); 
		colnames(M) = Rcpp::wrap( vInfo->ID );
	    rownames(M) = Rcpp::wrap( vInfo->sampleNames );  

		chunk = DataChunk<Rcpp::NumericMatrix, VariantInfo>( M, *vInfo );

		return ret;
	}


	virtual bool getNextChunk( DataChunk<vector<double>, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		chunk = DataChunk<vector<double>, VariantInfo>( matDosage, *vInfo );

		return ret;
	}




	private:
	View::UniquePtr view = nullptr; 
	size_t number_of_variants;
	size_t number_of_samples = 0;
	vector<string> sampleNames;
	map<size_t, size_t> requestedSamplesByIndexInDataIndex;
	VariantInfo *vInfo = nullptr;
	vector<double> probs;	
	vector<double> matDosage;
	size_t max_entries_per_sample = 3;		


	bool getNextChunk_helper(){	

		Dimension data_dimension = Dimension( number_of_variants, number_of_samples, max_entries_per_sample ) ;
		Dimension ploidy_dimension = Dimension( number_of_variants, number_of_samples ) ;

		IntegerVector ploidy = IntegerVector( ploidy_dimension, NA_INTEGER ) ;
		LogicalVector phased = LogicalVector( number_of_variants, NA_LOGICAL ) ;

		string SNPID, rsid, chromosome;
		genfile::bgen::uint32_t position;
		vector<string> alleles;

		// Iterate through variants
		for( size_t j = 0; j < number_of_variants; j++ ) {

			// read variant information
			view->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ;
		
			// store variant info
			vInfo->addVariant(chromosome, position, rsid, alleles[0], alleles[1] );

			// read genotype probabilities into DataSetter object
			DataSetter setter(
				&ploidy, ploidy_dimension,
				&probs, data_dimension,
				&phased,
				j,
				requestedSamplesByIndexInDataIndex
			);
			
			view->read_genotype_data_block( setter ) ;		
		}

		// Convert to dosage values stored in vector<double>
		//------------------

		// use probs to create Cube 
		cube C(probs.data(), number_of_variants, number_of_samples, max_entries_per_sample, true, true);
	
		vec dsg = {0,1,2}; // weight alleles by dosage

		// Using Tmp matrix
		// mat M(number_of_samples, number_of_variants);
		// for(int j=0; j<number_of_variants; j++){
		// 	M.col(j) = C.row_as_mat(j).t() * dsg;
		// }
	    // chunk = DataChunk<arma::mat, VariantInfo>( M, *vInfo );

		// compute dosage from Cube
		// copy results of each variant to vector<double>
		for(int j=0; j<number_of_variants; j++){
			vec v = C.row_as_mat(j).t() * dsg;
			memcpy(matDosage.data() + number_of_samples*j, v.memptr(), number_of_samples*sizeof(double));
		}

		return false;
	}
};


}

#endif