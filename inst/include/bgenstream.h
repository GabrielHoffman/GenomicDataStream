/***********************************************************************
 * @file		bgenstream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		reads a BGEN into matrix in chunks, storing variants in columns
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef BGEN_STREAM_H_
#define BGEN_STREAM_H_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef USE_EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 

#include <string>

#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"

#include <boost/algorithm/string.hpp>

#include <VariantInfo.h>
#include <GenomicDataStream.h>
#include <GenomicRanges.h>
#include "load.h"

using namespace std;
using namespace arma;
using namespace genfile::bgen;


namespace GenomicDataStreamLib {

/** Construct view of BGEN file using index to subset variants
 * based on region or variant id
 * @param filename path to BGEN file
 * @param index_filename path to index for BGEN file
 * @param gr `GenomicRanges` of intervals
 * @param rsids vector<string> of variant ids
 */ 
genfile::bgen::View::UniquePtr construct_view(
	const string & filename,
	const string & index_filename,
	const GenomicRanges & gr,
	const vector<string> & rsids = vector<string>()) {

	using namespace genfile::bgen ;

	View::UniquePtr view = View::create( filename ) ;

	if( gr.size() > 0){		
		IndexQuery::UniquePtr query = IndexQuery::create( index_filename ) ;
		for( int i = 0; i < gr.size(); i++ ) {
			query->include_range( IndexQuery::GenomicRange( gr.get_chrom(i) , gr.get_start(i), gr.get_end(i) ) ) ;
		}
		query->include_rsids( rsids ) ;
		query->initialise() ;
		view->set_query( query ) ;
	}
	return view ;
}

/** bgenstream reads a BGEN into an matrix in chunks, storing variants in columns.  Applies filtering for specified samples and genome region. 
 * 
*/
class bgenstream : 
	public GenomicDataStream {
	public:

	/** constructor
	*/
	bgenstream(const Param & param) : GenomicDataStream(param) {

		// Initialize genomic regions
		GenomicRanges gr( param.regions );

		// Initialize view and filter variants
		view = construct_view( param.file, param.file + ".bgi", gr ) ;

		// number of variants after filtering 
		n_variants_total = view->number_of_variants() ;

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

	/** Get number of columns in data matrix
	 */ 
	int n_samples(){
		return number_of_samples;
	}

	virtual bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		arma::mat M(matDosage.data(), number_of_samples, vInfo->size(), false, true);

	    chunk = DataChunk<arma::mat, VariantInfo>( M, *vInfo );

		return ret;
	}

	#ifdef USE_EIGEN
	virtual bool getNextChunk( DataChunk<Eigen::MatrixXd, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(matDosage.data(), number_of_samples, vInfo->size());

		chunk = DataChunk<Eigen::MatrixXd, VariantInfo>( M, *vInfo );

		return ret;
	}
	#endif

	#ifdef USE_RCPP
	virtual bool getNextChunk( DataChunk<Rcpp::NumericMatrix, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Rcpp::NumericMatrix M(number_of_samples, vInfo->size(), matDosage.data()); 
		colnames(M) = Rcpp::wrap( vInfo->ID );
	    rownames(M) = Rcpp::wrap( vInfo->sampleNames );  

		chunk = DataChunk<Rcpp::NumericMatrix, VariantInfo>( M, *vInfo );

		return ret;
	}
	#endif


	virtual bool getNextChunk( DataChunk<vector<double>, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		chunk = DataChunk<vector<double>, VariantInfo>( matDosage, *vInfo );

		return ret;
	}

	private:
	View::UniquePtr view = nullptr; 
	size_t number_of_samples = 0;
	vector<string> sampleNames;
	map<size_t, size_t> requestedSamplesByIndexInDataIndex;
	VariantInfo *vInfo = nullptr;
	vector<double> probs;	
	vector<double> matDosage;
	size_t max_entries_per_sample = 3;
	int n_variants_total;	
	int variant_idx_start = 0;

	bool getNextChunk_helper(){	

		// clear data, but keep allocated capacity
		matDosage.clear();
		vInfo->clear();

		// number of variants in this chunk
		size_t chunkSize = min(param.chunkSize, n_variants_total - variant_idx_start);

		// if no variants remain, return false
		if( chunkSize == 0) return false;

		vector<int> data_dimension;
		data_dimension.push_back(chunkSize);
		data_dimension.push_back(number_of_samples);
		data_dimension.push_back(max_entries_per_sample);

		vector<int> ploidy_dimension;
		ploidy_dimension.push_back( chunkSize ) ;
		ploidy_dimension.push_back( number_of_samples ) ;

		vector<int> ploidy(ploidy_dimension[0]*ploidy_dimension[1], numeric_limits<int>::quiet_NaN());
		vector<bool> phased( chunkSize, numeric_limits<bool>::quiet_NaN() ) ;

		string SNPID, rsid, chromosome;
		genfile::bgen::uint32_t position;
		vector<string> alleles;

		// Iterate through variants
		size_t k = 0;
		for( size_t j = variant_idx_start; j < variant_idx_start + chunkSize; j++ ) {

			// read variant information
			view->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ;
		
			// store variant info
			vInfo->addVariant(chromosome, position, rsid, alleles[0], alleles[1] );

			// read genotype probabilities into DataSetter object
			DataSetter setter(
				&ploidy, ploidy_dimension,
				&probs, data_dimension,
				&phased,
				k++,
				requestedSamplesByIndexInDataIndex
			);
			
			view->read_genotype_data_block( setter ) ;		
		}
		// increament starting position to beginning of next chunk
		variant_idx_start += chunkSize;

		// Convert to dosage values stored in vector<double>
		//------------------

		// use probs to create Cube 
		cube C(probs.data(), chunkSize, number_of_samples, max_entries_per_sample, true, true);
	
		vec dsg = {0,1,2}; // weight alleles by dosage

		// compute dosage from Cube
		// copy results of each variant to vector<double>
		for(int j=0; j<chunkSize; j++){
			// compute dosages
			vec v = C.row_as_mat(j).t() * dsg;

			// replace missing with mean
			if( param.missingToMean ) nanToMean( v );				

	    	// save vector in matDosage
			memcpy(matDosage.data() + number_of_samples*j, v.memptr(), number_of_samples*sizeof(double));
		}

		return true;
	}
};


// Using Tmp matrix
// mat M(number_of_samples, chunkSize);
// for(int j=0; j<chunkSize; j++){
// 	M.col(j) = C.row_as_mat(j).t() * dsg;
// }
// chunk = DataChunk<arma::mat, VariantInfo>( M, *vInfo );


}

#endif