/***********************************************************************
 * @file		bgenstream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		reads a BGEN into matrix in chunks, storing variants in columns
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef BGEN_STREAM_H_
#define BGEN_STREAM_H_

#ifndef DISABLE_EIGEN
#include <Eigen/Sparse>
#endif 

#include <string>

#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"

#include <boost/algorithm/string.hpp>

#include "VariantInfo.h"
#include "GenomicDataStream_virtual.h"
#include "GenomicRanges.h"
#include "bgen_load.h"

using namespace std;
using namespace arma;
using namespace genfile::bgen;

namespace gds {

/** Construct view of BGEN file using index to subset variants
 * based on region or variant id
 * @param filename path to BGEN file
 * @param index_filename path to index for BGEN file
 * @param gr `GenomicRanges` of intervals
 * @param rsids vector<string> of variant ids
 */ 
static genfile::bgen::View::UniquePtr construct_view(
	const string & filename,
	const string & index_filename,
	const GenomicRanges & gr,
	const vector<string> & rsids = vector<string>()) {

	if( ! filesystem::exists( filename ) ){
		throw runtime_error("File does not exist: " + filename);
	}
	if( ! filesystem::exists( index_filename ) ){
		throw runtime_error("File does not exist: " + index_filename);
	}

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

	bgenstream() {}

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
	View::UniquePtr view = nullptr; 
	size_t number_of_samples = 0;
	vector<string> sampleNames;
	map<size_t, size_t> requestedSamplesByIndexInDataIndex;
	VariantInfo *vInfo = nullptr;
	vector<double> probs;	
	vector<double> matDosage;
	size_t max_entries_per_sample = 4;
	int n_variants_total;	
	int variant_idx_start = 0;

	bool getNextChunk_helper(){	

		// clear data, but keep allocated capacity
		matDosage.clear();
		vInfo->clear();

		// number of variants in this chunk
		int chunkSize = min(param.chunkSize, n_variants_total - variant_idx_start);
		chunkSize = max(chunkSize, 0);

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
	
		// weight alleles by dosage
		// With max_entries_per_sample = 4, the unphased coding is 
		//  AA/AB/BB/NULL so use weights 0/1/2/0 to conver to dosage
		//  since the last entry doesn't encode valid information
		// When phased, the coding is [a1 a2] / [a1 a2]
		//  so use weights 0/1/0/1
		// vec w_unph = {0,1,2,0}; 
		// vec w_ph = {0,1,0,1}; 
		// compute dosages with weights depend on phasing
		// w = phased[j] ? w_ph : w_unph;
		// dsg = C.row_as_mat(j).t() * w;
		// BUT !!!
		// in unphased data, they last entry can be NaN
		//    and NaN * 0 is still NaN
		// so need to drop the last entry _manually_

		vec v, dsg;
		vec w_unph = {0,1,2}; 
		vec w_ph = {0,1,0,1}; 
		mat m;

		// compute dosage from Cube
		// copy results of each variant to vector<double>
		for(int j=0; j<chunkSize; j++){
			if( phased[j] ){
				dsg = C.row_as_mat(j).t() * w_ph;
			}else{		
				// extract columns 0,1,2
				// skip 3rd
				m = C.row_as_mat(j).t();		
				dsg = m.cols(0,2) * w_unph;
			}

			// replace missing with mean
			if( param.missingToMean ) nanToMean( dsg );				

	    	// save vector in matDosage
			memcpy(matDosage.data() + number_of_samples*j, dsg.memptr(), number_of_samples*sizeof(double));
		}

		return true;
	}
};

}

#endif
