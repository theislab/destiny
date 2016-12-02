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

IntegerVector rank(NumericVector x) {
	return match(x, clone(x).sort());
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> icor2_no_censor(
	const Eigen::SparseMatrix<double> dists,
	NumericMatrix imputed_data,
	const Function callback,
	bool use_rank = false
) {
	const int n = dists.rows();
	
	NumericMatrix data;
	if (!use_rank) {
		data = imputed_data;
	} else {
		data = NumericMatrix(imputed_data.nrow(), imputed_data.ncol());
		for (int i=0; i<data.nrow(); i++) {
			data(i, _) = rank(imputed_data(i, _));
		}
	}
	
	typedef Eigen::SparseMatrix<double> M;
	M d2(dists);
	
	for (int i=0; i<n; ++i) {
		for (M::InnerIterator it(d2, i); it; ++it) {
			const double d = 1 - cor(data(i, _), data(it.index(), _));
			it.valueRef() = d*d;
		}
		if (i%1000 == 0)
			callback(i+1);
	}
	callback(n);
	
	return d2;
}
