/***********************************************************************
 * @file		GenomicRanges.h.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Store genomic ranges of type chrom:start-end
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef GENOMIC_RANGES_H_
#define GENOMIC_RANGES_H_

#include <vector>
#include <string>
#include <stdlib.h>

#include "utils.h"

using namespace std;

namespace gds {

/** Store genomic ranges of type chrom:start-end
 */
class GenomicRanges {

	public:
	/** Constructor from vectors of chrom, start, end
	 */ 
	GenomicRanges( const vector<string> &chrom, const vector<uint32_t> &start, const vector<uint32_t> &end) : 
		chrom(chrom), start(start), end(end)
		{ }

	/** Constructor from string of chr:start-end
	 */
	GenomicRanges( const vector<string> &regions ){
		initialize( regions );
	}

	/** Constructor from string of delimited chr:start-end,chr:start-end for delim "\t,\n"
	 */
	GenomicRanges( const string &regionString ){

		vector<string> regions;
		boost::split(regions, regionString, boost::is_any_of("\t,\n"));
		initialize( regions );
	}

	/** Accessors
	 */
	const string get_chrom(const int &i) const { return chrom[i];}
	const uint32_t get_start(const int &i) const { return start[i];}
	const uint32_t get_end(const int &i) const { return end[i];}

	const int size() const { return end.size();}


	/** Evaluate if position is within [start,end] inclusive
	 */ 
	const bool isWithin( const uint32_t &start, const uint32_t &end, const uint32_t &position) const {

	    return (position >= start) && (position <= end); 
	}
	const bool isWithin( const string &chr, const uint32_t &position) const {

		bool found = false;
	    for (int i = 0; i < chrom.size(); i++) {
	        if (chrom[i].compare(chr) == 0 ) {
	        	if( isWithin( start[i], end[i], position) ){
	        		found = true;
	        		break;
	        	}
	        }
	    }

	    return found;
	}

	/** get indeces of entries in position that are found in Genomic ranges. Currently quadratic time
	 */
	const vector<int> getWithinIndeces( vector<string> chr,
										vector<uint32_t> position){

		vector<int> indeces;
		for(int i=0; i<chr.size(); i++){
			if( isWithin( chr[i], position[i]) ){
				indeces.push_back(i);
			}
		}

		return indeces;
	}

	const vector<int> getWithinIndeces( vector<string> chr,
										vector<string> position){

		auto &tmp = convert_to_uint32_t( position );
		return getWithinIndeces(chr, tmp);
	}


	private:
	vector<string> chrom;
	vector<uint32_t> start;
	vector<uint32_t> end;

	void initialize(vector<string> regions){

		// if only entry is "." or "", 
		// don't add any regions and return early
		if( regions.size() == 1 && regions[0].compare(".") == 0) return;
		if( regions.size() == 1 && regions[0].compare("") == 0) return;

		// Remove duplicate entries, but preserve element order. 
		removeDuplicates( regions );

		vector<string> reg, pos;
		long p0, p1;

		for(auto const & it : regions){

			// parse
			boost::split(reg, it, boost::is_any_of(":"));
			boost::split(pos, reg[1], boost::is_any_of("-"));

			// if only start give, set end to same value
			if( pos.size() == 1 ) pos.push_back( pos[0] );

			// convert strings to long int
			p0 = atol(pos[0].c_str());
			p1 = atol(pos[1].c_str());

			// check region
			if( (reg.size() != 2) || (pos.size() > 2) || (p1 < p0) || (p0<0) || (p1<0)){ 
				throw logic_error("Not a valid genomic range: " + it);
			}
			if( !isOnlyDigits(pos[0]) || !isOnlyDigits(pos[1])){
				throw logic_error("Not a valid genomic range: " + it);				
			}

			// save regions
			chrom.push_back( reg[0] );
			start.push_back( p0 );
			end.push_back( p1 );
		}	
	}

	/** Convert vector<string> to vector<uint32_t>
	 */ 
	const vector<uint32_t> convert_to_uint32_t( const vector<string> &v ){
		vector<uint32_t> output(0, v.size());

		for (auto &s : v) {
		    stringstream parser(s);
		    uint32_t x = 0;
		    parser >> x;
		    output.push_back(x);
		}
		return output;
	}
};


}

#endif