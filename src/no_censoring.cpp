// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include "utils.h"

using namespace Rcpp;


// [[Rcpp::export]]
Eigen::SparseMatrix<double> icor2_no_censor(
		const IntegerVector i,
		const IntegerVector j,
		const int n,
		const NumericMatrix imputed_data,
		const Function callback,
		const bool use_rank = false
) {
	NumericMatrix data;
	if (!use_rank) {
		data = clone(imputed_data);
	} else {
		data = NumericMatrix(imputed_data.nrow(), imputed_data.ncol());
		for (int r=0; r<data.nrow(); r++) {
			data(r, _) = rank(imputed_data(r, _));
		}
	}
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(i.size());
	
	for (int k=0; k<i.size(); k++) {
		const double d = 1 - cor(data(i[k], _), data(j[k], _));
		triplets.push_back(T(i[k], j[k], d*d));
		
		if (k%1000 == 0)
			callback(k+1);
	}
	callback(i.size());
	
	Eigen::SparseMatrix<double> d2(n, n);
	d2.setFromTriplets(triplets.begin(), triplets.end());
	return d2;
}
