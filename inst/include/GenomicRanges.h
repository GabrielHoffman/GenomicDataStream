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

namespace GenomicDataStreamLib {


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

	private:
	vector<string> chrom;
	vector<uint32_t> start;
	vector<uint32_t> end;

	void initialize(vector<string> regions){

		// if only entry is ".", don't add any regions and return early
		if( regions.size() == 1 && regions[0].compare(".") == 0) return;

		// Remove duplicate entries, but preserve element order. 
		removeDuplicates( regions );

		vector<string> reg, pos;

		for(auto const & it : regions){

			// parse
			boost::split(reg, it, boost::is_any_of(":"));
			boost::split(pos, reg[1], boost::is_any_of("-"));

			// if only start give, set end to same value
			if( pos.size() == 1 ) pos.push_back( pos[0] );

			// check region
			if( (reg.size() != 2) || (pos.size() > 2) || (pos[1] < pos[2])){ 
				throw logic_error("Not a valid genomic range: " + it);
			}
			if( !isOnlyDigits(pos[0]) || !isOnlyDigits(pos[1])){
				throw logic_error("Not a valid genomic range: " + it);				
			}

			// save regions
			chrom.push_back( reg[0] );
			start.push_back( stoul(pos[0]) );
			end.push_back( stoul(pos[1]) );
		}	
	}
};

}

#endif