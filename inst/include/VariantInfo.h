/***********************************************************************
 * @file        VariantInfo.h
 * @author      Gabriel Hoffman
 * @email       gabriel.hoffman@mssm.edu
 * @brief       Store variant information and sample names
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef VARIANT_INFO_H
#define VARIANT_INFO_H

#include <vector>
#include "DataInfo.h"
#include "utils.h"

using namespace std;

namespace gds {

/** Store variant information and sample names
 * 
 */
class VariantInfo : 
    public DataInfo {

    public:

    VariantInfo() {};

    /** constructor to intitialize sampleNames
     * @param sampleNames retained samples after filtering
     */ 
    VariantInfo( const vector<string> &sampleNames) : 
    	sampleNames(sampleNames) {};

    /** add information for one variant at a time
     */ 
    void addVariant(const string &chr,
                    const int & pos, 
                    const string &id, 
                    const string &allele1, 
                    const string &allele2){

        CHROM.push_back( chr );
        POS.push_back( pos );
        ID.push_back( id );
        A1.push_back( allele1 );
        A2.push_back( allele2 );
    }

    /** add information for many variant at a time
     */ 
    void addVariants(const vector<string> &chr,
                    const vector<string> &pos, 
                    const vector<string> &id, 
                    const vector<string> &allele1, 
                    const vector<string> &allele2){

        int len = chr.size();

        if( pos.size() != len || 
            id.size() != len ||
            allele1.size() != len || 
            allele1.size() != len){
            throw logic_error("All vectors must have same length");
        }

        CHROM.insert(CHROM.end(), 
            chr.begin(), 
            chr.end());

        // convert pos strings to ints
        vector<int> posInt;
        for(auto &x: pos){
          posInt.push_back( stoi(x) );
        }

        POS.insert(POS.end(), 
            posInt.begin(), 
            posInt.end());

        ID.insert(ID.end(),
            id.begin(),
            id.end());

        A1.insert(A1.end(),
            allele1.begin(),
            allele1.end());

        A2.insert(A2.end(),
            allele2.begin(),
            allele2.end());
    }

    /** append variants in a new VariantInfo to the end of the current object
     */ 
    void append( const VariantInfo & vInfo){

    	CHROM.insert(CHROM.end(), 
    		vInfo.CHROM.begin(), 
    		vInfo.CHROM.end());

    	POS.insert(POS.end(), 
    		vInfo.POS.begin(), 
    		vInfo.POS.end());

    	ID.insert(ID.end(),
            vInfo.ID.begin(),
            vInfo.ID.end());

    	A1.insert(A1.end(),
            vInfo.A1.begin(),
            vInfo.A1.end());

    	A2.insert(A2.end(),
            vInfo.A2.begin(),
            vInfo.A2.end());
    }

    /** Retain only variants with indeces stored in idx
     */ 
    void retainVariants( const vector<unsigned int> &idx){
        CHROM = subset_vector( CHROM, idx );
        POS = subset_vector( POS, idx );
        ID = subset_vector( ID, idx );
        A1 = subset_vector( A1, idx );
        A2 = subset_vector( A2, idx );
    }

    /** Clear vectors storing variant information, but leave sampleNames
     */ 
    void clear(){
        CHROM.clear();
        POS.clear();
        ID.clear();
        A1.clear();
        A2.clear();
    }

    /** get genome interval as a string
     */ 
    string getInterval(){

      if( POS.size() == 0) return "";

      // start
      string reg = CHROM[0] + ":" + to_string(POS[0]);
  
      // if there is a change in chromosome
      for(int i=1; i<POS.size(); i++){
        if( CHROM[i] != CHROM[i-1]){
          // end interval and start new one
          reg += "-" + to_string(POS[i-1]);
          reg += "," + CHROM[i] + ":" + to_string(POS[i]);
        }
      }

      reg += "-" + to_string(POS[POS.size() -1]);

      return reg;
    }

    // private:
    vector<string> sampleNames;
    vector<string> CHROM;
    vector<int> POS;
    vector<string> A1, A2;
};

}

#endif
