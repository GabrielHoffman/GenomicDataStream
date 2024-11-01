
// #include <RcppArmadillo.h>


// #include "Rtatami.h"
// #include <vector>
// #include <algorithm>


// // Not necessary in a package context, it's only used for this vignette:
// // [[Rcpp::depends(beachmat, assorthead)]]

// // [[Rcpp::export(rng=false)]]
// Rcpp::NumericVector column_sums(Rcpp::RObject initmat) {
//     Rtatami::BoundNumericPointer parsed(initmat);
//     const auto& ptr = parsed->ptr;

//     auto NR = ptr->nrow();
//     auto NC = ptr->ncol();
//     std::vector<double> buffer(NR);
//     Rcpp::NumericVector output(NC);
//     auto wrk = ptr->dense_column();

//     for (int i = 0; i < NC; ++i) {
//         auto extracted = wrk->fetch(i, buffer.data());
//         output[i] = std::accumulate(extracted, extracted + NR, 0.0);
//     }

//     return output;
// }

// //   sourceCpp('test.cpp')

// std::vector<double> output;

// // [[Rcpp::export(rng=false)]]
// arma::mat get_rows(Rcpp::RObject initmat, const std::vector<int> &ridx, const int &nthreads = 1) {
//     Rtatami::BoundNumericPointer parsed(initmat);
//     const auto& ptr = parsed->ptr;

//     auto NR = ridx.size();
//     auto NC = ptr->ncol();
//     output.reserve(NC*NR);

//     std::vector<double> buffer(NC);
//     auto wrk = ptr->dense_row();
//     for (int i = 0; i < ridx.size(); i++) {
//         auto extracted = wrk->fetch(ridx[i], buffer.data());
//         std::memcpy(output.data() + NC*i, extracted, NC*sizeof(double));
//     }

//     // map std::vector to arma::mat
//     arma::mat M(output.data(), NC, NR, false, true);

//     return M;
// }


