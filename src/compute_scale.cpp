#include "Rtatami.h"
#include "tatami_stats/tatami_stats.hpp"
#include <cmath>
#include <algorithm>

// fix issue that file and API changes variance -> variances
// fixed: https://chatgpt.com/c/6aa30b0c-fb5c-83e9-9e4d-e32632144184
// tatami_stats_compat.hpp

#ifndef MY_PACKAGE_TATAMI_STATS_COMPAT_HPP
#define MY_PACKAGE_TATAMI_STATS_COMPAT_HPP

#if __has_include(<tatami_stats/variances.hpp>)

#include <tatami_stats/variances.hpp>

// old API
namespace my_tatami_variance = tatami_stats::variances;

#elif __has_include(<tatami_stats/variance.hpp>)

#include <tatami_stats/variance.hpp>

// new API -- adjust this alias to the namespace used by the version
// you are targeting.
namespace my_tatami_variance = tatami_stats::variance;

#else

#error "Unsupported version of tatami_stats"

#endif

#endif


// [[Rcpp::export(rng=false)]]
Rcpp::List compute_center_and_scale(Rcpp::RObject mat, int nthreads) {
  Rtatami::BoundNumericPointer bound(mat);
  const auto& ptr = bound->ptr;
  auto NR = ptr->nrow();
  auto NC = ptr->ncol();

  Rcpp::NumericVector center(NC), scale(NC);
  double* cptr = static_cast<double*>(center.begin());
  double* sptr = static_cast<double*>(scale.begin());

  // Handling edge cases.
  if (NR <= 1) {
    if (NR == 0) {
      std::fill(center.begin(), center.end(), R_NaReal);
    } else {
      if (ptr->prefer_rows()) {
        auto iptr = ptr->dense_row()->fetch(0, cptr);
        tatami::copy_n(iptr, NC, cptr);
      } else {
        auto iptr = ptr->dense_column()->fetch(0, cptr);
        tatami::copy_n(iptr, NR, cptr);
      }
    }
    std::fill(scale.begin(), scale.end(), R_NaReal);
    return Rcpp::List::create(
      Rcpp::Named("center") = center, 
      Rcpp::Named("scale") = scale
    );
  }

  if (ptr->prefer_rows()) {
    if (ptr->sparse()) {
      tatami::parallelize([&](size_t, int start, int len) -> void {
        auto ext = tatami::consecutive_extractor<true>(ptr.get(), true, 0, NR, start, len);
        std::vector<double> vbuffer(len);
        std::vector<int> ibuffer(len);

        std::vector<double> tmp_means(len), tmp_vars(len);
        my_tatami_variance::RunningSparse<double, double, int> runner(len, tmp_means.data(), tmp_vars.data(), false, start);
        for (int r = 0; r < NR; ++r) {
          auto range = ext->fetch(r, vbuffer.data(), ibuffer.data());
          runner.add(range.value, range.index, range.number);
        }
        runner.finish();

        std::copy(tmp_means.begin(), tmp_means.end(), cptr + start);
        for (auto& v : tmp_vars) {
          v = std::sqrt(v);
        }
        std::copy(tmp_vars.begin(), tmp_vars.end(), sptr + start);
      }, NC, nthreads);

    } else {
      tatami::parallelize([&](size_t, int start, int len) -> void {
        auto ext = tatami::consecutive_extractor<false>(ptr.get(), true, 0, NR, start, len);
        std::vector<double> buffer(len);

        std::vector<double> tmp_means(len), tmp_vars(len);
        my_tatami_variance::RunningDense<double, double, int> runner(len, tmp_means.data(), tmp_vars.data(), false);
        for (int r = 0; r < NR; ++r) {
          auto ptr = ext->fetch(r, buffer.data());
          runner.add(ptr);
        }
        runner.finish();

        std::copy(tmp_means.begin(), tmp_means.end(), cptr + start);
        for (auto& v : tmp_vars) {
          v = std::sqrt(v);
        }
        std::copy(tmp_vars.begin(), tmp_vars.end(), sptr + start);
      }, NC, nthreads);
    }

  } else {
    if (ptr->sparse()) {
      tatami::parallelize([&](size_t, int start, int len) -> void {
        tatami::Options opt;
        opt.sparse_extract_index = false;
        auto ext = tatami::consecutive_extractor<true>(ptr.get(), false, start, len, opt);
        std::vector<double> vbuffer(NR);
        for (int c = start, end = start + len; c < end; ++c) {
          auto range = ext->fetch(c, vbuffer.data(), NULL);
          auto paired = my_tatami_variance::direct(range.value, range.number, NR, false);
          cptr[c] = paired.first;
          sptr[c] = std::sqrt(paired.second);
        }
      }, NC, nthreads);

    } else {
      tatami::parallelize([&](size_t, int start, int len) -> void {
        auto ext = tatami::consecutive_extractor<false>(ptr.get(), false, start, len);
        std::vector<double> buffer(NR);
        for (int c = start, end = start + len; c < end; ++c) {
          auto ptr = ext->fetch(c, buffer.data());
          auto paired = my_tatami_variance::direct(ptr, NR, false);
          cptr[c] = paired.first;
          sptr[c] = std::sqrt(paired.second);
        }
      }, NC, nthreads);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("center") = center, 
    Rcpp::Named("scale") = scale
  );
}
