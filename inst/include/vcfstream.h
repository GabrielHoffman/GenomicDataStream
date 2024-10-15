/***********************************************************************
 * @file		vcfstream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		vcfstream reads a VCF into a matrix in chunks, storing variants in columns
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef VCF_STREAM_H_
#define VCF_STREAM_H_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef USE_EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif 



#include <string>
#include <vcfpp.h>

#include "VariantInfo.h"
#include "GenomicDataStream.h"
#include "utils.h"

using namespace std;
using namespace vcfpp;

namespace GenomicDataStreamLib {

/** vcfstream reads a VCF/BCF into an matrix in chunks, storing variants in columns.  Applies filtering for specified samples and genome region. 
 * 
*/
class vcfstream : 
	public GenomicDataStream {
	public:

	/** constructor initilizing with parameter values
	*/
	vcfstream(const Param & param) : GenomicDataStream(param) {
	
 		// initialize
		reader = new BcfReader( param.file );

		validRegions.reserve(param.regions.size());
		validRegions.clear();

		// check status of each region
		// retain only valid, non-empty regions in validRegions
		for(const string& region : param.regions){

			switch( reader->getStatus( region ) ){
				case 1: // region is vaild and not empty
				validRegions.push_back(region);	
				break;

				case 0: // the region is valid but empty.
				break;

				case -1: // there is no index file found.
				throw runtime_error("Could not retrieve index file");
				break;

				case -2: // the region is not valid
				throw runtime_error("region was not found: " + region );
				break;
			}
		}

		// initialize iterator
		itReg = validRegions.begin();

		reader->setRegion( *itReg );
		reader->setSamples( param.samples );

		// Initialize record with info in header
		record = new BcfRecord( reader->header ); 

		// 1: int; 2: float; 3: string; 0: error;
		fieldType = reader->header.getFormatType(param.field);

		if( fieldType == 3 && param.field.compare("GT"))
		throw std::runtime_error("field GT is the only supported string type");				

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
	~vcfstream(){
		if( reader != nullptr) delete reader;
		if( record != nullptr) delete record;
		if( vInfo != nullptr) delete vInfo;
	}

	/** Get number of columns in data matrix
	 */ 
	int n_samples(){
		return reader->nsamples;
	}

	bool getNextChunk( DataChunk<arma::mat, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		// mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem = true, strict = false)
		bool copy_aux_mem = false; // create read-only matrix without re-allocating memory
		arma::mat M(matDosage.data(), reader->nsamples, vInfo->size(), copy_aux_mem, true);

	    chunk = DataChunk<arma::mat, VariantInfo>( M, *vInfo );

		return ret;
	}

	#ifdef USE_EIGEN
	bool getNextChunk( DataChunk<Eigen::MatrixXd, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Eigen::MatrixXd M = Eigen::Map<Eigen::MatrixXd>(matDosage.data(), reader->nsamples, vInfo->size());

		chunk = DataChunk<Eigen::MatrixXd, VariantInfo>( M, *vInfo );

		return ret;
	}
	#endif

	#ifdef USE_RCPP
	bool getNextChunk( DataChunk<Rcpp::NumericMatrix, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		Rcpp::NumericMatrix M(reader->nsamples, vInfo->size(), matDosage.data()); 
		colnames(M) = Rcpp::wrap( vInfo->ID );
	    rownames(M) = Rcpp::wrap( vInfo->sampleNames );  

		chunk = DataChunk<Rcpp::NumericMatrix, VariantInfo>( M, *vInfo );

		return ret;
	}
	#endif

	bool getNextChunk( DataChunk<vector<double>, VariantInfo> & chunk){

		// Update matDosage and vInfo for the chunk
		bool ret = getNextChunk_helper();

		chunk = DataChunk<vector<double>, VariantInfo>( matDosage, *vInfo );

		return ret;
	}

	/** Concatenate variant identifiers 
	 * 
	 * @param record storing current variant
	 * */
	static string variantToString( const BcfRecord &record ) {
		string s = record.CHROM() + ":" + to_string(record.POS()) + " " + record.ID() + " " + record.REF() + " " + record.ALT();
		return s;
	}

	private:
	BcfReader *reader = nullptr;
	BcfRecord *record = nullptr;
	VariantInfo *vInfo = nullptr;
	vector<string>::iterator itReg;
	vector<string> validRegions;

	bool continueIterating = true;
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
	// set using reserve() to set initial capacity so avoid re-alloc
	vector<double> matDosage;	

	bool getNextChunk_helper(){

		// if end of file reached, return false
		if( ! continueIterating ) return continueIterating;

		// clear data, but keep allocated capacity
		matDosage.clear();
		vInfo->clear();

		// loop thru variant, updating the count each time
		unsigned int j;
		for(j=0; j < param.chunkSize; j++){

			// get next variant
			// if false, reached end of region
			if( ! reader->getNextVariant( *record ) ){

				// else go to next region
				itReg++;

				// if this was the last region
				// set continueIterating so false is retured at next call to 
				// getNextChunk_helper()
				// then break since no data left
				if( itReg == validRegions.end()){
					continueIterating = false;
					break;
				}

				// else
				// initialize the record for this region
				reader->setRegion( *itReg );
				reader->getNextVariant( *record ); 
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

					// get GT as int's with vector that is twice as long
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
		return true;
	}
};


} // end namespace

#endif