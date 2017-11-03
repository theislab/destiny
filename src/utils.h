#ifndef _UTILS_H
#define _UTILS_H

#include <Rcpp.h>

using namespace Rcpp;

double cor(const NumericVector v1, const NumericVector v2);
IntegerVector rank(const NumericVector x);
NumericMatrix rank_mat(const NumericMatrix x);

#endif
