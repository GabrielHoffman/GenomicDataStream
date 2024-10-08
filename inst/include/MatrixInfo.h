/***********************************************************************
 * @file        MatrixInfo.h
 * @author      Gabriel Hoffman
 * @email       gabriel.hoffman@mssm.edu
 * @brief       Store matrix information
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef MATRIX_INFO_H
#define MATRIX_INFO_H

#include <vector>

using namespace std;

namespace GenomicDataStreamLib {


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


/** For Rcpp::NumericMatrix with MatrixInfo, set the row and col names
 */ 
static void setRowColNames( Rcpp::NumericMatrix &M, const MatrixInfo &info){
    Rcpp::rownames(M) = info.get_rownames();
    Rcpp::colnames(M) = info.get_colnames();
}



/** For all other datatypes, do nothing
 */ 
template<typename matType, typename infoType >
static void setRowColNames( matType &M, const infoType &info){
} 



}

#endif