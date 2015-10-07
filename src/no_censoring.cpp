// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

double cor(const NumericVector v1,
           const NumericVector v2) {
	
	if (v1.size() != v2.size()) stop("v1 needs to be of same size as v2");
	
	const int n = v1.size();
	
	double v1sum, v2sum, v12sum, v1sqr_sum, v2sqr_sum;
	v1sum = v2sum = v12sum = v1sqr_sum = v1sqr_sum = 0;
	
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

void validate_n_k(const NumericMatrix nn_dist, const int n, const int k) {
	char buffer[100];
	const char* error = "number of %s of nn_index (%d) and nn_dist (%d) must be equal";
	if (n != nn_dist.nrow()) { sprintf(buffer, error, "rows", n, nn_dist.nrow()); stop(buffer); }
	if (k != nn_dist.ncol()) { sprintf(buffer, error, "cols", k, nn_dist.ncol()); stop(buffer); }
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> d2_no_censor(
	const IntegerMatrix nn_index,
	const NumericMatrix nn_dist,
	const Function callback
) {
	const int n = nn_index.nrow();
	const int k = nn_index.ncol();
	
	validate_n_k(nn_dist, n, k);
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(k * n);
		
	for (int i=0; i<n; i++) {
		for (int j=0; j<k; j++) {
			const int    nn_i = nn_index(i, j) - 1;
			const double d    = nn_dist (i, j);
			triplets.push_back(T(i, nn_i, d*d));
		}
		if (i%1000 == 0)
			callback(i+1);
	}
	callback(n);
	
	Eigen::SparseMatrix<double> d2(n, n);
	d2.setFromTriplets(triplets.begin(), triplets.end());
  return d2;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> icos2_no_censor(
	const IntegerMatrix nn_index,
	NumericMatrix imputed_data,
	const Function callback
) {
	const int n = nn_index.nrow();
	const int k = nn_index.ncol();
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(k * n);
	
	for (int i=0; i<n; i++) {
		for (int j=0; j<k; j++) {
			const int nn_i = nn_index(i, j) - 1;
			const double d = 1 - cor(imputed_data(i, _), imputed_data(nn_i, _));
			triplets.push_back(T(i, nn_i, d*d));
		}
		if (i%1000 == 0)
			callback(i+1);
	}
	callback(n);
	
	Eigen::SparseMatrix<double> d2(n, n);
	d2.setFromTriplets(triplets.begin(), triplets.end());
	return d2;
}
