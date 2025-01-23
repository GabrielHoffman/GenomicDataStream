/***********************************************************************
 * @file		VariantSet.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Set of chrom/position/ID searchable with genome interval
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef VARIANT_SET_H
#define VARIANT_SET_H

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>  

#include "GenomicRanges.h"

using namespace std;

namespace gds {

/** Store genomic position of variant and its index 
 */ 
struct point {
	point(const size_t & position, const size_t & index) :
		position(position), index(index) {}

	point(const size_t & position) :
		position(position) {}

	size_t position;
	size_t index;
};

/** Overload less-than operator for use by sorting functions is_sorted(), lower_bound(), upper_bound()
 */ 
static bool operator<( const point &a, const point &b){
	return a.position < b.position;
}

/** Store position and ID, sorted by position.  Get indeces of variants with a query interval using binary search in O(log(N)) time.
 */
class VariantSet {

	public:
	VariantSet( const vector<string> &chrom, const vector<size_t> &position){

		// for each variant
		// insert (position[i], i) into chromosome hash
		for(int i=0; i<chrom.size(); i++){
			map[chrom[i]].push_back( point(position[i], i) );
		}

		// for each chrom
		// check that positions are sorted
		for( auto &v: map){
			if( ! is_sorted(v.second.begin(), v.second.end()) ){
				throw logic_error("Positions are not sorted in " + v.first);
			}
		}
	}

	/** Get indeces of variants within the query interval using binary search 
	*/
	vector<size_t> getIndeces( const string &chrom, const size_t &start, const size_t &end){

		auto vBegin = map[chrom].begin();
		auto vEnd = map[chrom].end();

		// get iterator to first element in the interval
		auto it1 = lower_bound(vBegin, vEnd, point(start)); 

		// get iterator to last element in the interval
		auto it2 = upper_bound(vBegin, vEnd, point(end));  

		vector<size_t> indeces;

		// walk from lower to upper bound
		// saving index at each step
		while( it1 != it2){
			indeces.push_back( it1->index );
			it1++;
		}

		return indeces;
	}

	/** Get indeces of variants within set of query intervals.  Use  binary search for each interval
	*/
	vector<size_t> getIndeces( const GenomicRanges &gr ){

		vector<size_t> indeces;

		// for each genome interval
		for(int i=0; i<gr.size(); i++){
			// get indeces of variants within this interval
			vector<size_t> idx = getIndeces(gr.get_chrom(i), gr.get_start(i), gr.get_end(i));

			// insert into vector
   			indeces.insert(indeces.end(), idx.begin(), idx.end());
		}

		return indeces;
	}

	private:
	unordered_map<string, vector<point> > map;
};

}

#endif
