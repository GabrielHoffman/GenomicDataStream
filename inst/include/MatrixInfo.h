/***********************************************************************
 * @file        MatrixInfo.h
 * @author      Gabriel Hoffman
 * @email       gabriel.hoffman@mssm.edu
 * @brief       Store matrix information
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef MATRIX_INFO_H
#define MATRIX_INFO_H

#ifndef DISABLE_DELAYED_STREAM

#include <vector>
#include <algorithm>
#include <Rcpp.h>
#include "DataInfo.h"

using namespace std;

namespace gds {

class MatrixInfo :
    public DataInfo {
    public:
    MatrixInfo() {}
    MatrixInfo(const vector<string> &ID) :
        DataInfo(ID) {}

    void setRowNames( const vector<string> &rn, const int &start, const int &end){

        ID.clear();
        for(int i=start; i<end; i++){
            ID.push_back(rn[i]);
        }
    }
};




}

#endif
#endif
