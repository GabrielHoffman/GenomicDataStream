

#ifdef ARMA
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif


using namespace arma;


#include <vector>
#include <random>
#include <algorithm>


#ifndef UTILS_H_
#define UTILS_H_

namespace GenomicDataStreamLib {

/** Compute dosage values from vector of GT stored as int.  Sum adjacent values to get dosage
* @param v for n samples, vector of length 2*n where dosage is computed as v[2*i] + v[2*i+1];
* @param missingToMean if true, set missing values to the mean dosage value.  if false, set to NaN
*/
static vector<double> intToDosage( const vector<int> &v, const bool &missingToMean) {

    // store and return result
    vector<double> res( v.size() / 2.0);
    vector<int> missing;

    // initialize
    int runningSum = 0, nValid = 0;
    double value;

    // for each entry in result
    // use two adjacent values
    for(int i=0; i<res.size(); i++){
        value = v[2*i] + v[2*i+1];

        // -9 is the missing value, so -18 is diploid
        if( value == -18){ 
            // if missing, set to NaN
            value = std::numeric_limits<double>::quiet_NaN();
            missing.push_back(i);
        }else{
            // for computing mean
            runningSum += value;
            nValid++;
        }

        // set dosage value
        res[i] = value;
    }

    // mean excluding NaNs
    double mu = runningSum / (double) nValid;

    // if missing values should be set to mean
    if( missingToMean ){
        // for each entry with a missing value, set to mean
        for(const int& i : missing) res[i] = mu;
    }

    return res;
}

/** Remove duplicate entries, but preserve element order.  Based on https://stackoverflow.com/questions/12200486/how-to-remove-duplicates-from-unsorted-stdvector-while-keeping-the-original-or
* Uses C++11
 */ 
template<typename T>
static size_t removeDuplicates(vector<T>& vec){
    unordered_set<T> seen;

    auto newEnd = remove_if(vec.begin(), vec.end(), [&seen](const T& value)
    {
        if (seen.find(value) != end(seen))
            return true;

        seen.insert(value);
        return false;
    });

    vec.erase(newEnd, vec.end());

    return vec.size();
}



/** Compute sum of each column
 * @param X matrix
 */ 
static arma::vec colSums( const arma::mat &X){

    // row vector of 1's
    arma::rowvec ONE(X.n_rows, fill::ones);

    // matrix multiplication to get sums
    arma::mat tmp = ONE * X;

    return arma::conv_to<arma::vec>::from( tmp );
}


/** Center and scale columns of matrix
 * @param X matrix with features as columns
 * @param center center columns
 * @param center scale columns by sd
 * 
 */
static void standardize( arma::mat &X, const bool &center = true, const bool &scale = true ){
    
    double sqrt_rdf = sqrt(X.n_rows - 1.0);

    // if center, subtract mean of each column
    // if scale, divide by sd of each column
    // Note, norm() does not center the column
    //   this give results consistent with base::scale()
    //   when scale is FALSE
    for(size_t j=0; j<X.n_cols; j++){
        if( center ) X.col(j) -= mean(X.col(j));
        if( scale )  X.col(j) /= norm(X.col(j)) / sqrt_rdf;
    }
}

}


#endif