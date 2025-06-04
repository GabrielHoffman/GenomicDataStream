#include "Rtatami.h"
#include "tatami_stats/tatami_stats.hpp"
#include <cmath>
#include <algorithm>


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
                tatami_stats::variances::RunningSparse<double, double, int> runner(len, tmp_means.data(), tmp_vars.data(), false, start);
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
                tatami_stats::variances::RunningDense<double, double, int> runner(len, tmp_means.data(), tmp_vars.data(), false);
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
                    auto paired = tatami_stats::variances::direct(range.value, range.number, NR, false);
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
                    auto paired = tatami_stats::variances::direct(ptr, NR, false);
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
