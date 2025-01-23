/***********************************************************************
 * @file		pgenstream.h
 * @author		Gabriel Hoffman
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

#include "VariantInfo.h"
#include "GenomicDataStream_virtual.h"
#include "GenomicRanges.h"
#include "DataTable.h"
#include "pgen/RPgenReader.h"
#include "VariantSet.h"

using namespace std;
using namespace arma;


/*      TODO

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

/** pgenstream reads a PGEN into an matrix in chunks, storing variants in columns.  Applies filtering for specified samples and genome region. 
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

		// Parse PSAM/FAM file of sample identifiers
		// evaluate subsetting of samples
    process_samples();

		// Read data file (pgen/bed)
		pg = new RPgenReader();
  	pg->Load(param.file, pvar, n_samples_psam, sampleIdx1);

		// Initialize vector with capacity to store nVariants
		// Note, this allocates memory but does not change .size()
		// After j variants have been inserted, only entries up to j*nsamples are populated
		//  the rest of the vector is allocated doesn't have valid data
		int n = 1e6 * param.initCapacity / (double) (sizeof(double) * number_of_samples);
		matDosage.reserve( n );
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
	size_t number_of_samples = 0;
	int n_requested_variants = 0;
	int currentIdx = 0;
	vector<double> matDosage;
	vector<size_t> varIdx;
	VariantInfo *vInfo = nullptr;
	RPgenReader *pg = nullptr;
	RPvar *pvar = nullptr;
	DataTable dt;
	// unordered_map linking ID to index
	unordered_map<string,int> map_dt_id;
	string fileIdx;
	int n_samples_psam;
    vector<int> sampleIdx1;
	FileType genoFileType;

	bool getNextChunk_helper2 (){

		// clear data, but keep allocated capacity
		matDosage.clear();
		vInfo->clear();

		// number of variants in this chunk
		int chunkSize = min(param.chunkSize, n_requested_variants - currentIdx);
		chunkSize = max(chunkSize, 0);

		// if no variants remain, return false
		if( chunkSize == 0) return false;

		// indeces of variants in chunk
		vector<int> varIdx_sub = 
						{varIdx.begin() + currentIdx, 
						 varIdx.begin() + currentIdx + chunkSize};

		// read dosage into matDosage using 					
		// 1-based indeces
		vector<int> varIdx_sub1(varIdx_sub);
		for(int &i : varIdx_sub1) i++; 

		pg->ReadList( matDosage, varIdx_sub1, param.missingToMean);

		// populate vInfo
  		string id, a1, a2,chrom = "";
  		int pos = -1;

  		// if PGEN
		if( genoFileType == PGEN){

	  		// Get variant info from pvar 
	  		for(auto i: varIdx_sub){
	  			id = pvar->GetVariantId(i);
	  			a1 = pvar->GetAlleleCode(i,0);
	  			a2 = pvar->GetAlleleCode(i,1);

	  			// find chrom, pos given id
			    int idx = map_dt_id[id];
			    chrom = dt["CHROM"][idx];
			    pos = atoi(dt["POS"][idx].c_str());

	  			vInfo->addVariant(chrom, pos, id, a1, a2);
	  		}
	  	}else{
	  		// Get variant info from DataTable from BIM
	  		for(auto i: varIdx_sub){
	  			chrom = dt["CHROM"][i];
	  			pos = atoi(dt["POS"][i].c_str());
	  			id = dt["ID"][i];
	  			a1 = dt["REF"][i];
	  			a2 = dt["ALT"][i];

	  			vInfo->addVariant(chrom, pos, id, a1, a2);
	  		}
	  	}

  		// increment current index
  		currentIdx += varIdx_sub.size();

  		return true;
	}

	bool getNextChunk_helper(){return getNextChunk_helper2();}

	/*  Parse PVAR file of variant positions and IDs
	*/
	void process_variants(){

		// if file is PGEN
    	if( genoFileType == PGEN ){
			// Name of .pvar file based on replacing .pgen$
			fileIdx = regex_replace(param.file, regex("pgen$"), "pvar");

			// Read .pvar file into DataTable
			// column names are define by line starting with "#CHROM"
			// lines before this are ignored
			dt = DataTable( fileIdx, "#CHROM" );

			genoFileType = PGEN;

		// if file is BED
		}else if( genoFileType == PBED ){

			// Name of .bim file based on replacing .pgen$
			fileIdx = regex_replace(param.file, regex("bed$"), "bim");

			// Read BIM file with no headerKey
			dt = DataTable( fileIdx );
			dt.setColNames({"CHROM", "ID", "CM", "POS", "ALT", "REF"});

			genoFileType = PBED;

		}else{
			throw logic_error("Not valid genotype file extension: " + param.file);
		}

		// populate unordered_map linking ID to index
		// to allow fast search
		for(int i=0; i<dt["ID"].size(); i++){
			map_dt_id.emplace(dt["ID"][i], i); 
		}

		// Initialize genomic regions
		// from delimited string
		GenomicRanges gr( param.regions );

		// if not empty
		if( gr.size() != 0){
			// get indeces of entries in .pvar located 
			// within param.regions
			// Search is linear time for each interval
			// varIdx = gr.getWithinIndeces( dt["CHROM"], cast_elements<size_t>(dt["POS"]) );

			// Search is log time for each interval
			VariantSet vs(dt["CHROM"], cast_elements<size_t>(dt["POS"]));
			varIdx = vs.getIndeces( gr );

		}else{
			// else
			// set entries to seq(0, dt.nrows()-1)
			varIdx.resize(dt.nrows());
			iota(begin(varIdx), end(varIdx), 0); 
		}

		// total number of requested variants
		n_requested_variants = varIdx.size();
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
			dt2 = DataTable( fileSamples );

			// set column names
			vector<string> names = {"FID", "IID", "PID", "MID", "SEX", "ALT", "PHENO"};
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
		}
	}

};

}

#endif
