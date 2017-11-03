#include "utils.h"

using namespace Rcpp;

double cor(const NumericVector v1,
           const NumericVector v2) {
	
	const int n = v1.size();
	if (n != v2.size()) stop("v1 needs to be of same size as v2");
	
	double v1sum, v2sum, v12sum, v1sqr_sum, v2sqr_sum;
	v1sum = v2sum = v12sum = v1sqr_sum = v2sqr_sum = 0;
	
	for (int i=0; i<n; i++) {
		v1sum += v1[i];
		v2sum += v2[i];
		v12sum += v1[i] * v2[i];
		v1sqr_sum += v1[i] * v1[i];
		v2sqr_sum += v2[i] * v2[i];
	}
	
	const double num = n*v12sum - v1sum*v2sum;
	const double deno = (n*v1sqr_sum - v1sum*v1sum) * (n*v2sqr_sum - v2sum*v2sum);
	
	return num / sqrt(deno);
}


IntegerVector rank(const NumericVector x) {
	return match(x, clone(x).sort());
}

// [[Rcpp::export]]
NumericMatrix rank_mat(const NumericMatrix x) {
	NumericMatrix ranked = NumericMatrix(x.nrow(),  x.ncol());
	for (int r=0; r<ranked.nrow(); r++) {
		ranked(r, _) = rank(x(r, _));
	}
	return ranked;
}
