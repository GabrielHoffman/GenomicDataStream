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

#include "VariantInfo.h"
#include "GenomicDataStream_virtual.h"
#include "GenomicRanges.h"
#include "DataTable.h"
#include "pgen/RPgenReader.h"

using namespace std;
using namespace arma;


/*      TODO

1) custom psam file
2) plink 1 support
3) empty headerKey
4) GenomicRanges search is quadratic, doesn't use sorting
5) subset samples
6) fill in chrom and pos info

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
		
		// Name of .pvar file based on replacing .gen$
		string fileIdx;
		fileIdx = regex_replace(param.file, regex("pgen$"), "pvar");
		if( ! filesystem::exists( fileIdx ) ){
			fileIdx = regex_replace(param.file, regex("pgen$"), "pvar.gz");
		}

		// Rcpp::Rcout << "Load fileIdx to DataTable..." << endl;
		dt = DataTable( fileIdx, "#CHROM" );

		// Initialize genomic regions
		GenomicRanges gr( param.regions );

		if( gr.size() != 0){
			// get indeces of entries in param.fileIdx located 
			// with in param.regions
			varIdx = gr.getWithinIndeces( dt["CHROM"], dt["POS"]);
		}else{
			// set entries to seq(0, dt.nrows()-1)
			varIdx.resize(dt.nrows());
			iota(begin(varIdx), end(varIdx), 0); 
		}

		n_variants_total = varIdx.size();

		// Rcpp::Rcout << "Variants retained: " << n_variants_total << endl;

		// Read index file (pvar/bim)
		// Rcpp::Rcout << "Load fileIdx..." << endl;
		pvar = new RPvar();
        pvar->Load(fileIdx);
       	// Rcpp::Rcout << "# variants: " << pvar->GetVariantCt()<< endl;

  		// initialize varInfo with sample names
  		string fileSamples = regex_replace(param.file, regex("pgen$"), "psam");
		DataTable dt2(fileSamples, "#IID");
  		vector<string> SamplesNames = dt2["IID"];

		// Rcpp::Rcout << "Print SamplesNames..." << endl;
		// for(auto i: SamplesNames) Rcpp::Rcout << i << endl;

		vInfo = new VariantInfo( SamplesNames );
		number_of_samples = SamplesNames.size();

        vector<int> sample_subset_1based;
		int raw_sample_ct = SamplesNames.size();

		// Read data file (pgen/bed)
		// Rcpp::Rcout << "Load PGEN..." << endl;
		pg = new RPgenReader();
  		pg->Load(param.file, pvar, raw_sample_ct, sample_subset_1based);

		// Initialize vector with capacity to store nVariants
		// Note, this allocates memory but does not change .size()
		// After j variants have been inserted, only entries up to j*nsamples are populated
		//  the rest of the vector is allocated doesn't have valid data
		int n = 1e6 * param.initCapacity / (double) (sizeof(double) * SamplesNames.size());
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
	int n_variants_total = 0;
	int currentIdx = 0;
	vector<double> matDosage;
	vector<int> varIdx;
	VariantInfo *vInfo = nullptr;
	RPgenReader *pg = nullptr;
	RPvar *pvar = nullptr;
	DataTable dt;

	bool getNextChunk_helper2 (){

		// clear data, but keep allocated capacity
		matDosage.clear();
		vInfo->clear();

		// number of variants in this chunk
		int chunkSize = min(param.chunkSize, n_variants_total - currentIdx);
		chunkSize = max(chunkSize, 0);

		// if no variants remain, return false
		if( chunkSize == 0) return false;

		// indeces of variants in chunk
		vector<int> varIdx_sub = {varIdx.begin() + currentIdx, 
								  varIdx.begin() + currentIdx + chunkSize};

		// read dosage into matDosage using 					
		// 1-based indeces
		vector<int> varIdx_sub1(varIdx_sub);
		for(int &i : varIdx_sub1) i++; 
		pg->ReadList( matDosage, varIdx_sub1, param.missingToMean);

		// populate vInfo
  		string id, a1, a2;
  		string chrom = "";
  		int pos = -1;
  		for(auto i: varIdx_sub){
  			id = pvar->GetVariantId(i);
  			a1 = pvar->GetAlleleCode(i,0);
  			a2 = pvar->GetAlleleCode(i,1);

  			// find chrom, pos given id
  			auto it = find(dt["ID"].begin(), dt["ID"].end(), id);
		    int idx = distance(dt["ID"].begin(), it);
		    chrom = dt["CHROM"][idx];
		    pos = atoi(dt["POS"][idx].c_str());

  			vInfo->addVariant(chrom, pos, id, a1, a2);
  		}

  		// increment current index
  		currentIdx += varIdx_sub.size();

  		return true;
	}

	bool getNextChunk_helper(){return getNextChunk_helper2();}
};

}

#endif