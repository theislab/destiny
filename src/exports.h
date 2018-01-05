#ifndef _DESTINY_EXPORTS_H
#define _DESTINY_EXPORTS_H

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

Eigen::SparseMatrix<double> censoring_impl(const NumericMatrix data, const NumericVector sigmas, const Eigen::SparseMatrix<double> dists, const SEXP thr_or_null, const SEXP uncertain_or_null, const SEXP missing_or_null, const Function callback);
NumericMatrix predict_censoring_impl(const NumericMatrix data, const NumericMatrix data2, const double thr, const NumericVector uncertain, const NumericVector missing, const double sigma);
List knn_cross(const NumericMatrix data, const NumericMatrix query, const size_t k, const std::string distance);
List knn_asym(const NumericMatrix data, const size_t k, const std::string distance);
NumericMatrix rank_mat(const NumericMatrix x);

#endif
