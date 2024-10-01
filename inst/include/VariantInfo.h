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

using namespace std;


class MatrixInfo {

    public:
    MatrixInfo() {}

    MatrixInfo(const Rcpp::CharacterVector &rownames, const Rcpp::CharacterVector &colnames) :
        rownames(rownames), colnames(colnames) {}

    MatrixInfo( const Rcpp::RObject &mat){
        Rcpp::List lst;
        if( mat.hasAttribute("dimnames") ){
            lst = mat.attr("dimnames");
        }else if( mat.hasAttribute("Dimnames") ){
            lst = mat.attr("Dimnames");
        }
        rownames = lst[0];
        colnames = lst[1];
    }

    /** Accessor
     */
    Rcpp::CharacterVector get_rownames() const { 
        return rownames;
    }

    /** Accessor
    */
    Rcpp::CharacterVector get_colnames() const{ 
        return colnames;
    }

    private:
    Rcpp::CharacterVector rownames, colnames;
};



/** Store variant information and sample names
 * 
 */
class VariantInfo {

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
                    const int64_t & pos, 
                    const string &id, 
                    const string &ref, 
                    const string &alt){

        CHROM.push_back( chr );
        POS.push_back( pos );
        ID.push_back( id );
        REF.push_back( ref );
        ALT.push_back( alt );
    }

    /** number of variants
     */ 
    int size() const {
        return CHROM.size();
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

    	REF.insert(REF.end(),
                        vInfo.REF.begin(),
                        vInfo.REF.end());

    	ALT.insert(ALT.end(),
                        vInfo.ALT.begin(),
                        vInfo.ALT.end());
    }

    /** Clear vectors storing variant information, but leave sampleNames
     */ 
    void clear(){
        CHROM.clear();
        POS.clear();
        ID.clear();
        REF.clear();
        ALT.clear();
    }

    /** Add information about statistical test for each variant.
     * NOT COMPLETE
     * @param testResults stores results from analysis with fastlmmBatchDesign
     */ 
    // void cbind( const vector<fastlmm_result> & testResults){
        
    //     // place-holder
    // }

    // private:
    vector<string> sampleNames;
    vector<string> CHROM;
    vector<int64_t> POS;
    vector<string> ID, REF, ALT;
};



#endif