/***********************************************************************
 * @file		vcfstream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		vcfstream reads a VCF into an arma::mat in chunks, storing variants in columns
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef VCF_STREAM_H_
#define VCF_STREAM_H_

#include <RcppArmadillo.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <string>
#include <vcfpp.h>

#include <VariantInfo.h>
#include "/Users/gabrielhoffman/workspace/repos/GenomicDataStream/inst/include/GenomicDataStream.h"

using namespace std;
using namespace vcfpp;

/** vcfstream reads a VCF into an arma::mat in chunks, storing variants in columns.  Applies filtering for specified samples and genome region. 
 * 
*/
class vcfstream : 
	public GenomicDataStream {
	public:

	/** constructor
	 * @param file .vcf, .vcf.gz or .bcf file with tabix index
	 * @param field 'GT' for genotype strings, 'DS' for dosage, or another other field stored as an integer or float.  'GT' is the only string type supported
	 * @param region target in the format 'chr2:1-12345'. Setting region to "" includes all variants
	 * @param samples string of comma separated sample IDs to extract: "ID1,ID2,ID3"
	 * @param chunkSize number of variants to return per chunk
	 * @param missingToMean if true, set missing values to the mean dosage value.  if false, set to NaN
	 * @param initCapacity initial capacity of temporary vector to avoid re-alloc on insert.  Size is in Mb.
	*/
	vcfstream(const Param & param) : GenomicDataStream(param) {
	
 		// initialize
		reader = new BcfReader( param.file );

		// check if region is valid
		switch( reader->getStatus( param.region) ){
			
			case 1: // region is vaild and not empty

				reader->setRegion(param.region);
				reader->setSamples(param.samples);

				// Initialize record with info in header
				record = new BcfRecord( reader->header ); 

				// 1: int; 2: float; 3: string; 0: error;
				fieldType = reader->header.getFormatType(param.field);

				if( fieldType == 3 && param.field.compare("GT"))
					throw std::runtime_error("field GT is the only supported string type");				
				break;
			
			case 0: // the region is valid but empty.
 				break;
 			
			case -1: // there is no index file found.
				throw std::runtime_error("Could not retrieve index file");
				break;
			
			case -2: // the region is not valid
				throw std::runtime_error("region was not found! make sure the region format is correct");
				break;
		}

		// initialize varInfo with sample names
		vInfo = new VariantInfo( reader->SamplesName );

		// Initialize vector with capacity to store nVariants
		// Note, this allocates memory but does not change .size()
		// After j variants have been inserted, only entries up to j*nsamples are populated
		//  the rest of the vector is allocated doesn't have valid data
		int n = 1e6 * param.initCapacity / (double) (sizeof(double) * reader->nsamples);
		matDosage.reserve( n );
	}

	/** destructor
	 */ 
	virtual ~vcfstream(){
		if( reader != nullptr) delete reader;
		if( record != nullptr) delete record;
		if( vInfo != nullptr) delete vInfo;
	}


	virtual bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		// mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem = true, strict = false)
		bool copy_aux_mem = false; // create read-only matrix without re-allocating memory
		arma::mat M(matDosage.data(), reader->nsamples, vInfo->size(), copy_aux_mem, true);

		chunk = DataChunk<arma::mat, VariantInfo>( M, *vInfo );

		return ret;
	}

	virtual bool getNextChunk( DataChunk<Eigen::Map<Eigen::MatrixXd>, 
		VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Eigen::Map<Eigen::MatrixXd> M(matDosage.data(), reader->nsamples, vInfo->size());

		chunk = DataChunk<Eigen::Map<Eigen::MatrixXd>, VariantInfo>( M, *vInfo );

		return ret;
	}


	/** Concatenate variant identifiers 
	 * 
	 * @param record storing current variant
	 * */
	static string variantToString( const BcfRecord &record ){
		string s = record.CHROM() + ":" + to_string(record.POS()) + " " + record.ID() + " " + record.REF() + " " + record.ALT();
		return s;
	}


	/** Compute dosage values from vector of GT stored as int.  Sum adjacent values to get dosage
	* @param v for n samples, vector of length 2*n where dosage is computed as v[2*i] + v[2*i+1];
	* @param missingToMean if true, set missing values to the mean dosage value.  if false, set to NaN
	*/
	static vector<double> intToDosage( const vector<int> &v, const bool &missingToMean){

		// store and return result
		vector<double> res( v.size() / 2.0);
		vector<int> missing;

		// initialize
		int runningSum = 0, nValid = 0;
		double value;

		// for each entry in result
		// use two adjacent values
		for(int i=0; i<res.size(); i++){
			value = v[2*i] + v[2*i+1];

			// -9 is the missing value, so -18 is diploid
			if( value == -18){ 
				// if missing, set to NaN
				value = std::numeric_limits<double>::quiet_NaN();
				missing.push_back(i);
			}else{
				// for computing mean
				runningSum += value;
				nValid++;
			}

			// set dosage value
			res[i] = value;
		}

		// mean excluding NaNs
		double mu = runningSum / (double) nValid;

		// if missing values should be set to mean
		if( missingToMean ){
			// for each entry with a missing value, set to mean
			for(const int& i : missing) res[i] = mu;
		}

		return res;
	}

	private:
	BcfReader *reader = nullptr;
	BcfRecord *record = nullptr;
	VariantInfo *vInfo = nullptr;

	bool endOfFile = false;
	int fieldType;

	// store genotype values
	// type used based on fieldType
	// 1: int; 
	// 2: float; 
	// 3: string; 
	// 0: error;
	vector<int> values_int; 
	vector<float> values_fl; 

	vector<double> tmp;	

	// stores genotype dosage as doubles, with the next marker inserted at the end
	// NOTE that when current size is exceeded, .insert() reallocates memory
	// this can be slow 
	// set using reserve() to set initial capacity so avaoid re-alloc
	vector<double> matDosage;	



	bool getNextChunk_helper(){

		// clear data, but keep allocated capacity
		matDosage.clear();
		vInfo->clear();

		// loop thru variant, updating the count each time
		unsigned int j;
		for(j=0; j < param.chunkSize; j++){

			// get next variant
			// if false, set endOfFile and break
			if( ! reader->getNextVariant( *record ) ){
				endOfFile = true;
				break;
			}

			// populate genotype with the values of the current variant
			// If string, convert to dosage
			// use values vector based on fieldType
			switch(fieldType){
				case 1: // int					
					record->getFORMAT( param.field, values_int);
					matDosage.insert(matDosage.end(), values_int.begin(), values_int.end());
					break;
				case 2: // float				
					record->getFORMAT( param.field, values_fl);
					matDosage.insert(matDosage.end(), values_fl.begin(), values_fl.end());
					break;
				case 3: // string. Convert GT to doubles

					// check if site is multi-allelic
					if( record->isMultiAllelics() ){
						throw std::runtime_error("GT is not supported for multi-allelic site\n    " + variantToString(*record));
					}
					// check if site ploidy > 2
					if( record->ploidy() > 2 ){
						throw std::runtime_error("GT is not supported for site with ploidy > 2\n    " + variantToString(*record));
					}					

					// get GT as int's with a vector that is twice as long 				
					record->getGenotypes(values_int);
					tmp = intToDosage( values_int, param.missingToMean );
					matDosage.insert(matDosage.end(), tmp.begin(), tmp.end());
					break;
			}

			// store variant information 
			vInfo->addVariant(	record->CHROM(), 
								record->POS(), 
								record->ID(), 
								record->REF(), 
								record->ALT() );
		}

		// After j variants have been inserted, only entries up to j*nsamples are populated
		//  the rest of the vector is allocated doesn't have valid data.
		// Therefore, use entires based on .size(), not .capacity()
		return ! endOfFile;
	}
};




#endif