/***********************************************************************
 * @file        DataTable.h
 * @author      Gabriel Hoffman
 * @email       gabriel.hoffman@mssm.edu
 * @brief       Load delimited data table from file, and access using column names
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef DATA_TABLE_H
#define DATA_TABLE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <stdlib.h> 

using namespace std;

namespace gds {

class DataTable {

	public:

	DataTable(){}

	DataTable(const string &file, const string &headerKey, const char delim = '\t'){

		if( ! filesystem::exists( file ) ){
			throw logic_error("File does not exist: " + file);
		}

    	ifstream strm( file );
    	string line, header, value;

    	// Loop through lines until header start key is found
    	bool startHeader = false;

    	// for each line
		while (getline(strm, line)) {

			// for each column
		    stringstream ss(line);
		    while (getline(ss, header, delim)) {

		    	// if not started yet, and found start yet
		    	// set startHeader to true
		    	if( ! startHeader ){
		    		if( header.compare(headerKey) == 0 ){
		       			startHeader = true;
		    			// remove leading # from header key
		    			header = regex_replace(header, regex("^#"), "");
		       		}
		    	}

		    	// if start key already found, 
		    	// push header column
		    	if( startHeader ){
            		colNames.push_back(header);
		       	}
		    }

		    // if end of line where start key is found
		    // break
		    if( startHeader ) break;
		}

		// initialize data with column for each header entry
		for(int i=0; i<colNames.size(); i++){
			data.push_back( vector<string>() );
		}

		// Read data rows
		int lineIdx = 0;
		// for each row
		while (getline(strm, line)) {
		    stringstream ss(line);
		    int colIndex = 0;

		    // for each column
		    while (getline(ss, value, delim)) {
		    	if( colIndex > data.size() ) break;
		        data[colIndex++].push_back( value );
		    }

		    if( colIndex != colNames.size()){
				throw logic_error("Line " + to_string(lineIdx) + " is not valid");
		    }

		    lineIdx++;
		}

		strm.close();
	}

	~DataTable(){}

	int ncols(){
		return data.size();
	}	

	int nrows(){
		return data[0].size();
	}	

	vector<string> getColNames(){
		return colNames;
	}	

	const vector<string> getCol(const string &key) const {
		// search for key in colNames
		auto it = std::find(colNames.begin(), colNames.end(), key);

		vector<string> ret;
 
 		// if found
		if (it != colNames.end()) {
			// get index key was found at
	        int index = distance(colNames.begin(), it);
	        ret = data[index];
	    }else{
			throw logic_error("Column not found: " + key);
	    } 

	    return ret;
	}

	const vector<string> operator[](const string &key) const {
		return getCol( key );
	}
	
	protected:
	vector<string> colNames; 
	vector<vector<string> > data;
};

}

#endif