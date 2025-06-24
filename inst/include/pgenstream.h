/***********************************************************************
 * @file		pgenstream.h
 * @author	Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		reads a plink2/PGEN into matrix in chunks, storing variants in columns
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef PGEN_STREAM_H_
#define PGEN_STREAM_H_

#ifndef DISABLE_EIGEN
#include <Eigen/Sparse>
#endif 

#include <string>
#include <regex>
#include <unordered_map>
#include <algorithm>

#include "VariantInfo.h"
#include "GenomicDataStream_virtual.h"
#include "GenomicRanges.h"
#include "DataTable.h"
#include "pgen/RPgenReader.h"
#include "VariantSet.h"
#include "utils.h"

using namespace std;
using namespace arma;


/*			TODO

done 1) custom psam file
done 2) plink 1 support
done 3) empty headerKey
done 4) GenomicRanges search is quadratic, doesn't use sorting
	test VariantSet
	how does distance() work?
done 5) subset samples
done 6) fill in chrom and pos info
		use hash instead

SamplesNames from file, and get raw_sample_ct
- getNextChunk_helper() logic for chunk size
- PVar must read positions
- intersect with BED file
- define VarIdx from BED file


*/

namespace gds {

/** pgenstream reads a PGEN into an matrix in chunks, storing variants in columns.	Applies filtering for specified samples and genome region. 
 * 
*/
class pgenstream : 
	public GenomicDataStream {
	public:

	pgenstream() {}

	/** constructor
	*/
	pgenstream(const Param & param) : GenomicDataStream(param) {
		
		genoFileType = getFileType(param.file);

		// Parse PVAR/BIM file of variant positions and IDs
		// evaluate subsetting of variants
		process_variants();

		if( genoFileType == PGEN ){
			// Read index file (pvar)
			pvar = new RPvar();
			pvar->Load(fileIdx, true, true);
		}

		// if PBED, set missingToMean to TRUE
		missingToMean = (genoFileType == PBED) ? 
											true : param.missingToMean;

		// Parse PSAM/FAM file of sample identifiers
		// evaluate subsetting of samples
		process_samples();

		// Read data file (pgen/bed)
		pg = new RPgenReader();
		pg->Load(param.file, pvar, n_samples_psam, sampleIdx1);

		// Initialize vector with capacity to store nVariants
		// Note, this allocates memory but does not change .size()
		// After j variants have been inserted, only entries up to j*nsamples are populated
		//	the rest of the vector is allocated doesn't have valid data
		matDosage.reserve( n_samples() * param.chunkSize );		
	}

	/** destructor
	 */ 
	~pgenstream(){
		if( vInfo != nullptr) delete vInfo;
		if( pg != nullptr){
			pg->Close();
			delete pg;
		}
		if( pvar != nullptr){
			pvar->Close();
			delete pvar;
		}
	}

	/** setter
	 */
	void setRegions(const vector<string> &regions ) override {

		// Initialize genomic regions
		// from delimited string
		GenomicRanges gr( regions );
		// if not empty
		if( gr.size() != 0){

			// Search is log time for each interval
			varIdx = vs.getIndeces( gr );

		}else{
			varIdx.clear();
			for(int i=0; i<dt.nrows(); i++){
				varIdx.push_back(i);
			}
		}

		// total number of requested variants
		n_requested_variants = varIdx.size();

		// set current position in varIdx
		currentIdx = 0;
	}

	/** Get number of columns in data matrix
	 */ 
	int n_samples() override {
		return number_of_samples;
	}

	/** Get vector of sample names in order that the genotypes are extracted
	 */ 
	vector<string> getSampleNames() override {
		return vInfo->sampleNames;
	}

	/** get FileType of param.file
	 */ 
	string getStreamType() override {
		return toString( param.fileType);
	}

	bool getNextChunk( DataChunk<arma::mat> & chunk, const bool &useFilter = true) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();	

		if( useFilter ){
			// modifies matDosage and vInfo directly
			applyVariantFilter(matDosage, vInfo, number_of_samples, getMAF(), getMinVariance() );
		}
		
		arma::mat M(matDosage.data(), number_of_samples, vInfo->size(), false, true);
		
		chunk = DataChunk<arma::mat>( M, vInfo );

		return ret;
	}

	bool getNextChunk( DataChunk<arma::sp_mat> & chunk, const bool &useFilter = true) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		if( useFilter ){
			// modifies matDosage and vInfo directly
			applyVariantFilter(matDosage, vInfo, number_of_samples, getMAF(), getMinVariance() );
		}

		arma::mat M(matDosage.data(), number_of_samples, vInfo->size(), false, true);

		chunk = DataChunk<arma::sp_mat>( arma::sp_mat(M), vInfo );

		return ret;
	}

	#ifndef DISABLE_EIGEN
	bool getNextChunk( DataChunk<Eigen::MatrixXd> & chunk, const bool &useFilter = true) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		if( useFilter ){
			// modifies matDosage and vInfo directly
			applyVariantFilter(matDosage, vInfo, number_of_samples, getMAF(), getMinVariance() );
		}		

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(matDosage.data(), number_of_samples, vInfo->size());

		chunk = DataChunk<Eigen::MatrixXd>( M, vInfo );

		return ret;
	}

	bool getNextChunk( DataChunk<Eigen::SparseMatrix<double> > & chunk, const bool &useFilter = true) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		if( useFilter ){
			// modifies matDosage and vInfo directly
			applyVariantFilter(matDosage, vInfo, number_of_samples, getMAF(), getMinVariance() );
		}	

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(matDosage.data(), number_of_samples, vInfo->size());

		chunk = DataChunk<Eigen::SparseMatrix<double>>( M.sparseView(), vInfo );

		return ret;
	}
	#endif

	#ifndef DISABLE_RCPP
	bool getNextChunk( DataChunk<Rcpp::NumericMatrix> & chunk, const bool &useFilter = true) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		if( useFilter ){
			// modifies matDosage and vInfo directly
			applyVariantFilter(matDosage, vInfo, number_of_samples, getMAF(), getMinVariance() );
		}	

		Rcpp::NumericMatrix M(number_of_samples, vInfo->size(), matDosage.data()); 
		colnames(M) = Rcpp::wrap( vInfo->getFeatureNames() );
		rownames(M) = Rcpp::wrap( vInfo->sampleNames );	

		chunk = DataChunk<Rcpp::NumericMatrix>( M, vInfo );

		return ret;
	}
	#endif

	bool getNextChunk( DataChunk<vector<double> > & chunk, const bool &useFilter = true) override {

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		if( useFilter ){
			// modifies matDosage and vInfo directly
			applyVariantFilter(matDosage, vInfo, number_of_samples, getMAF(), getMinVariance() );
		}	

		chunk = DataChunk<vector<double> >( matDosage, vInfo );

		return ret;
	}

	private:
	size_t number_of_samples = 0;
	int n_requested_variants = 0;
	int currentIdx;
	vector<double> matDosage;
	vector<int> varIdx;
	VariantInfo *vInfo = nullptr;
	RPgenReader *pg = nullptr;
	RPvar *pvar = nullptr;
	DataTable dt;
	VariantSet vs;
	string fileIdx;
	int n_samples_psam;
	vector<int> sampleIdx1;
	FileType genoFileType;
	bool missingToMean;

	bool getNextChunk_helper (){

		// clear data, but keep allocated capacity
		matDosage.clear();
		vInfo->clear();

		// number of variants in this chunk
		int chunkSize = min(param.chunkSize, n_requested_variants - currentIdx);
		chunkSize = max(chunkSize, 0);

		// if no variants remain, return false
		if( chunkSize == 0) return false;

		// indeces of variants in chunk
		auto end = min( varIdx.begin() + currentIdx + chunkSize,  varIdx.end());
		vector<int> varIdx_sub = {varIdx.begin() + currentIdx, end};

		// read dosage into matDosage using 					
		// 1-based indeces
		// vector<int> varIdx_sub1(varIdx_sub);
		vector<int> varIdx_sub1(varIdx_sub);
		for(int &i : varIdx_sub1) i++;

		pg->ReadList( matDosage, varIdx_sub1, missingToMean);

		// Populate vInfo from DataTable
		// Looking up column in DataTable is slow
		// so do it once and process vector
		vInfo->addVariants(	subset_vector(dt["CHROM"], varIdx_sub),
												subset_vector(dt["POS"], varIdx_sub), 
												subset_vector(dt["ID"], varIdx_sub), 
												subset_vector(dt["REF"], varIdx_sub), 
												subset_vector(dt["ALT"], varIdx_sub));

		// increment current index
		currentIdx += varIdx_sub.size();

		return true;
	}

	// bool getNextChunk_helper(){return getNextChunk_helper2();}

	/*	Parse PVAR file of variant positions and IDs
	*/
	void process_variants(){

		// if file is PGEN
		if( genoFileType == PGEN ){
			// Name of .pvar file based on replacing .pgen$
			fileIdx = regex_replace(param.file, regex("pgen$"), "pvar");

			// Read .pvar file into DataTable
			// column names are defined by line starting with "#CHROM"
			// lines before this are ignored
			dt = DataTable( fileIdx, "#CHROM" );

		// if file is BED
		}else if( genoFileType == PBED ){

			// Name of .bim file based on replacing .pgen$
			fileIdx = regex_replace(param.file, regex("bed$"), "bim");

			// Read BIM file with no headerKey

			dt = DataTable( fileIdx );
			dt.setColNames({"CHROM", "ID", "CM", "POS", "ALT", "REF"});
		}else{
			throw logic_error("Not valid genotype file extension: " + param.file);
		}
	
		// initialize variant set
		vs = VariantSet(dt["CHROM"], stoi_vec(dt["POS"]));

		// Set genomic regions regions
		setRegions( param.regions );
	}

	void process_samples(){

		DataTable dt2;

		string fileSamples;

		// Get path to samples file 
		// if custom samples file is given
		if( param.fileSamples.compare("") != 0){
			fileSamples = param.fileSamples;
		}else{
			if( genoFileType == PGEN ){
				// Name of .psam file based on replacing .pgen$
					fileSamples = regex_replace(param.file, regex("pgen$"), "psam");
			}else if( genoFileType == PBED ){
				// Name of .fam file based on replacing .pgen$
				fileSamples = regex_replace(param.file, regex("bed$"), "fam");
			}
		}

		// Read sample file depending on extension
		if( regex_search(fileSamples, regex("psam$")) ){
			// Read .psam file into DataTable
			// column names are define by line starting with "#IID"
			// lines before this are ignored
			dt2 = DataTable(fileSamples, "#IID");

		}else if( regex_search(fileSamples, regex("fam$")) ){

			// Read BIM file with no headerKey
			// space delimited entries
			dt2 = DataTable( fileSamples, "");

			// set column names
			vector<string> names = {"FID", "IID", "PID", "MID", "SEX", "PHENO"}; 
			vector<string> names_sub(names.begin(), names.begin() + dt.ncols());
			dt2.setColNames(names_sub);
		}else{
			throw logic_error("Not valid sample file extension: " + fileSamples);
		}

		// Sample names from PSAM file
		vector<string> SamplesNames = dt2["IID"];
		n_samples_psam = SamplesNames.size();

		// indeces of samples to include
		vector<int> sampleIdx;

		// Filter samples
		// If param.samples contains entries 
		if( param.samples.compare("-") != 0 ){

			// get sample ids from param.samples
			vector<string> requestedSamples;

			// split delmited string into vector
			boost::split(requestedSamples, param.samples, boost::is_any_of("\t,\n"));

			// Use unordered_map linking sample id to index
			// for fast searching
			// sn = sampleNames
			unordered_map<string,int> map_sn;
			for(int i=0; i<SamplesNames.size(); i++){
				map_sn.emplace(SamplesNames[i], i); 
			}

			// For each requested sample ID, 
			// get its index in the PSAM file
			for( string & name : requestedSamples){
				if( map_sn.count(name) == 0){
					throw logic_error("Sample id not found: " + name);
				}
				sampleIdx.push_back( map_sn[name] );
			}
			sort(sampleIdx.begin(), sampleIdx.end());

			// Get the requested sample ID's
			// sorted according to the PSAM file
			vector<string> requestedSamples_ordered;
			for(int i : sampleIdx){
				requestedSamples_ordered.push_back( SamplesNames[i] );	
			}

			// populate VariantInfo with requested sample IDs
			// in the order from the PSAM file
			vInfo = new VariantInfo( requestedSamples_ordered );
			number_of_samples = requestedSamples_ordered.size();

			// convert to 1-based indeces
			sampleIdx1.assign(sampleIdx.begin(), sampleIdx.end());
			for(int &i : sampleIdx1) i++; 
		}else{
			vInfo = new VariantInfo( SamplesNames );
			number_of_samples = SamplesNames.size();

			// set sampleIdx1 to be seq(1, number_of_samples)
			sampleIdx1.resize(number_of_samples);
			iota(begin(sampleIdx1), end(sampleIdx1), 1); 
		}
	}

}; 

}

#endif
