// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>

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

IntegerVector rank(NumericVector x) {
	return match(x, clone(x).sort());
}

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
