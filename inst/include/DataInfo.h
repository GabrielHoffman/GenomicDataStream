/***********************************************************************
 * @file        DataInfo.h
 * @author      Gabriel Hoffman
 * @email       gabriel.hoffman@mssm.edu
 * @brief       Meta-data (i.e. feature names) from dataset read by GenomimcDataStream
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef DATA_INFO_H
#define DATA_INFO_H

#include <vector>
#include <string>

using namespace std;

namespace gds {

class DataInfo {
    public:
    DataInfo() {}
    DataInfo(const vector<string> &ID) :
        ID(ID) {}

    vector<string> getFeatureNames() const {
        return ID;
    }

    string getFeatureName(const int &i) const {
        return ID[i];
    }

    /** number of entries
     */ 
    int size() const {
        return ID.size();
    }

    protected:
    vector<string> ID;
};

}

#endif
