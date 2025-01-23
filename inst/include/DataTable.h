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

	/** Read file into DataTable
		column names are define by line starting with \code{headerKey}
		lines before this are ignored.
		Columns of DataTable can then be accessed by name
		Only string values a supported, elements can be converted afterward

		@param file file path
		@param headerKey column names are define by line starting with \code{headerKey}
		@param delim single character delimiter
	*/
	DataTable(const string &file, const string &headerKey="", const char delim = '\t'){

		if( ! filesystem::exists( file ) ){
			throw logic_error("File does not exist: " + file);
		}

		// open file
	    ifstream strm( file );
		string line, header, value;

		// if there is NO headerKey
		if( headerKey.compare("") == 0){

			// read first line, count number of columns
			// create colNames
			ifstream strm_tmp( file );

			// get first line
			getline(strm_tmp, line);

			// for each column
			int ncols = 1;
		    stringstream ss(line);
		    while (getline(ss, header, delim)) {
				colNames.push_back("col" + to_string(ncols++));
			}
			strm_tmp.close();

		// if there is a headerKey
		}else{
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

			if( ! startHeader ){
				throw logic_error("Header key not found: " + headerKey);
			}
		}

		// initialize data with column for each header entry
		for(int i=0; i<colNames.size(); i++){
			data.push_back( vector<string>() );
		}

		// Read data rows
		// data is a vector of columns
		// each column is a vector of strings
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

	const int ncols() const {
		return data.size();
	}	

	const int nrows() const {
		return data[0].size();
	}	

	const vector<string> getColNames() const {
		return colNames;
	}	

	void setColNames(const vector<string> &names) {
		if( names.size() != colNames.size()){
			throw logic_error("setColNames: new and old names must have the same number of entries");
		}

  		colNames.assign( names.begin(), names.end());
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

	/** print DataTable to ostream
	 */ 
	void print(ostream& out, const string &delim = "\t") const {

		// print column names
		for(int i=0; i<colNames.size()-1; i++){
			out << colNames[i] << delim;
		}
		out << colNames[colNames.size()-1] << endl;

		// for each row
		for(int r=0; r<this->nrows(); r++){
			for(int i=0; i<data.size()-1; i++){
				out << data[i][r] << delim;
			}
			out << data[data.size()-1][r] << endl;
		}
	}
	
	protected:
	vector<string> colNames; 
	vector<vector<string> > data;
};

}

#endif
