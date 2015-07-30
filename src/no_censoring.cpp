// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::SparseMatrix<double> d2_no_censor(
	IntegerMatrix nn_index,
	NumericMatrix nn_dist,
  Function callback
) {
	int n = nn_index.nrow();
	int k = nn_index.ncol();
	
	char buffer[100];
	const char* error = "number of %s of nn_index (%d) and nn_dist (%d) must be equal";
	if (n != nn_dist.nrow()) { sprintf(buffer, error, "rows", n, nn_dist.nrow()); stop(buffer); }
	if (k != nn_dist.ncol()) { sprintf(buffer, error, "cols", k, nn_dist.ncol()); stop(buffer); }
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(k * n);
		
	for (int i=0; i<n; i++) {
		for (int j=0; j<k; j++) {
			int    nn_i = nn_index(i, j) - 1;
			double d    = nn_dist (i, j);
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
