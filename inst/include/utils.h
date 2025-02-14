
#include <armadillo>

#include <vector>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <regex>

#include <boost/algorithm/string.hpp>

#ifndef UTILS_H_
#define UTILS_H_

namespace gds {

/** Compute dosage values from vector of GP stored as double or float.  Sum adjacent values to get dosage
* @param v for n samples, vector of length 3*n where dosage is computed as c(v[3*i], v[3*i+1], v[3*i+2]) %*% c(0,1,2)
* @param missingToMean if true, set missing values to the mean dosage value.  if false, set to NaN
*/

template<typename T>
static vector<double> GP_to_dosage( const vector<T> &v, const bool &missingToMean) {
    vector<double> res( v.size() / 3.0);
    vector<int> missing;

    // initialize
    int runningSum = 0, nValid = 0;
    double value;

    // for each entry in result
    // use two adjacent values
    for(int i=0; i<res.size(); i++){
        // compute dosage from genotype probabilties
        value = v[3*i]*0 + v[3*i+1]*1 + v[3*i+2]*2;

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
    arma::rowvec ONE(X.n_rows, arma::fill::ones);

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


/** if string contains only digits, return true.  Else false
 */
static bool isOnlyDigits(const std::string& s){
    int n = count_if(s.begin(), s.end(),
                         [](unsigned char c){ return isdigit(c); } 
                        );
    return( n == s.size());
}


/** Replace nan values with mean
 */ 
static void nanToMean( arma::vec & v){
    // get indeces of finite elements
    arma::uvec idx = arma::find_finite(v);

    // if number of finite elements is less than the total
    if( idx.n_elem < v.n_elem ){
        // compute mean from finite elements
        double mu = arma::mean( v.elem(idx));

        // replace nan with mu
        v.replace(arma::datum::nan, mu);
    }
}


/** enum to indicate file type for genetics files
 */ 
typedef enum {
    VCF,
    VCFGZ,
    BCF,
    BGEN,
    PGEN,
    PBED,
    OTHER
} FileType;

/** return string from enum FileType
 */ 
static string toString( FileType x){

    switch(x){
        case VCF:   return "vcf";
        case VCFGZ:   return "vcf.gz";
        case BCF:   return "bcf";
        case BGEN:   return "bgen";
        case PGEN:   return "pgen";
        case PBED:   return "bed";
        case OTHER:   return "other";
        default:   return "other";
    }
}



/** Use regex to get file type from file name
 */ 
static FileType getFileType( const string &file ){

    FileType ft = OTHER;

    if( regex_search( file, regex("\\.vcf$")) ){
        ft = VCF;
    }else if( regex_search( file, regex("\\.vcf\\.gz$")) ){
        ft = VCFGZ;
    }else if( regex_search( file, regex("\\.bcf$")) ) {
        ft = BCF;
    }else if( regex_search( file, regex("\\.bgen$")) ){
        ft = BGEN;
    }else if( regex_search( file, regex("\\.pgen$")) ){
        ft = PGEN;
    } if( regex_search( file, regex("\\.bed$")) ){
        ft = PBED;
    }

    return ft;
}



/** Convert vector<string> to vector<T>
 */ 
template<typename T>
vector<T> cast_elements( const vector<string> &v ){
    vector<T> output(0, v.size());

    for (auto &s : v) {
        stringstream parser(s);
        T x = 0;
        parser >> x;
        output.push_back(x);
    }
    return output;
}


/** regionString is string of chr:start-end delim by "\t,\n"
 remove spaces, then split based on delim
 remove duplicate regions, but preserve order
 Note: regionString is copy by value, since boost::erase_all overwrites
 */
static vector<string> splitRegionString( string regionString){

    vector<string> regions;

    // regionString is string of chr:start-end delim by "\t,\n"
    // remove spaces, then split based on delim
    boost::erase_all(regionString, " ");
    boost::split(regions, regionString, boost::is_any_of("\t,\n"));

    // remove duplicate regions, but preserve order
    removeDuplicates( regions );

    return regions;
}


/** Return subset of x indicated by entries in ind: x[idx] 
*/
template<typename T>
static vector<T> subset_vector(const vector<T> &x, const vector<unsigned int> &idx){

    // initialize x_subset to have size idx.size()
    vector<T> x_subset;
    x_subset.reserve(idx.size());

    for(int i=0; i<idx.size(); i++){
        if( idx[i] > x.size() ){
            throw std::out_of_range("Index is out of bounds"); 
        }
        x_subset.push_back(x[idx[i]]);
    }

    return x_subset;
}

} // end namespace
#endif
