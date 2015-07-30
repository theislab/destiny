// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

#include <cmath>

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

double censor_pair(
	double c, double d,
	double sigma, double kt,
	double thr,
	double uncertain0, double uncertain1,
	double missing0, double missing1
) {
	bool use_d,
		one_uncertain  = (c == thr) || (d == thr),
		one_missing  = NumericVector::is_na(c) || NumericVector::is_na(d),
		one_of_each_invalid = one_uncertain && one_missing,
		both_valid = !one_uncertain && !one_missing;
	
	/*
	Rcout << "one uncertain=" << one_uncertain << " one missing=" << one_missing
				<< " one of each invalid=" << one_of_each_invalid << " both valid=" << both_valid << '\n';
	Rcout << "sigma=" << sigma << " kt=" << kt << '\n';
	Rcout << "before i=" << row_idx << " (" << c << "), j=" << nn_idx << " (" << d << "): " << x << '\n';
	*/
	
	//4 totally different models, depending on combinations
	// u: uncertain
	// m: missing
	// v: valid
	// 
	// |   | v | u | m |
	// |---|---|---|---|
	// | v | 1 | 2a| 2b|
	// | u | 2a| 4 | 3 |
	// | m | 2b| 3 | 4 |
	// 
	// 1. (valid×2): interference of two gaussians
	// 2. (invalid-valid): interference of box function and gaussian. a/b are parameters for the box
	// 3. (missing-uncertain): interference of two different box functions
	// 4. (invalid×2): two missing or two uncertain => 1
	
	if (both_valid) { //1
		return exp(-pow(c-d, 2) / (kt*2));
	} else { //at least one invalid
		if (one_of_each_invalid) { //3
			double uncertain_range = uncertain1 - uncertain0 + 2*sigma;
			return uncertain_range / (sqrt(uncertain_range) * sqrt(missing1 - missing0));
		} else if (one_uncertain || one_missing) { //2
			double m0, m1;
			if (one_uncertain && !one_missing) {
				m0  = uncertain0;
				m1  = uncertain1;
				use_d = c == thr;
			} else if (!one_uncertain && one_missing) {
				m0  = missing0;
				m1  = missing1;
				use_d = NumericVector::is_na(c);
			}
			
			double v = use_d ? d : c;
			
			return
					pow(M_PI*kt/2, -1./4)
				* sqrt(M_PI*kt/4)
				* ( std::erfc((m0-v) / sigma) - std::erfc((m1-v) / sigma) )
				/ sqrt(m1-m0);
			
			//Rcout << "m0=" << m0 << " m1=" << m1 << " use d=" << use_d << " x=" << x << " v=" << v << '\n';
		} else {// both are uncertain or both are missing => x *= 1 => noop
			return 1;
		}
	}
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> censoring_impl(
		NumericMatrix data,
		SEXP thr_or_null,
		SEXP uncertain_or_null,
		SEXP missing_or_null,
		double sigma,
		SEXP nns_or_null,
		Function callback
) {
	int
		n = data.nrow(),
		G = data.ncol();
	
	Rboolean
		no_nns = Rf_isNull(nns_or_null),
		no_threshold = Rf_isNull(thr_or_null),
		no_uncertain = Rf_isNull(uncertain_or_null),
		no_missing = Rf_isNull(missing_or_null);
	
	IntegerMatrix nns         = (no_nns)       ? IntegerMatrix(0, n) : IntegerMatrix(nns_or_null);
	NumericVector thr_v       = (no_threshold) ? NumericVector(0)    : NumericVector(thr_or_null);
	NumericMatrix uncertain_m = (no_uncertain) ? NumericMatrix(0, 2) : NumericMatrix(uncertain_or_null);
	NumericMatrix missing_m   = (no_missing)   ? NumericMatrix(0, 2) : NumericMatrix(missing_or_null);
	
	//the sparse matrix has about k·n elements (or k·(n+1)?)
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve((no_nns) ? (n*n) : (nns.ncol() * (n+1)));
	
	double kt = pow(sigma, 2);
	
	for (int g_idx=0; g_idx<uncertain_m.nrow(); g_idx++) {
		uncertain_m(g_idx, 0) -= sigma;
		uncertain_m(g_idx, 1) += sigma;
	}
	
	for (int g_idx=0; g_idx<missing_m.nrow(); g_idx++) {
		missing_m(g_idx, 0) -= sigma;
		missing_m(g_idx, 1) += sigma;
	}
	
	for (int row_idx=0; row_idx<n; row_idx++) {
		for (int j=0; j<nns.ncol(); j++) {
			int nn_idx;
			if (no_nns)
				nn_idx = j;
			else
				nn_idx = nns(row_idx, j) - 1;
			
			//Rcout << "i=" << row_idx << ", j=" << nn_idx << '\n';
			 
			double x = 1.;
		
			for (int g=0; g<G; g++) {
				double
					c = data(row_idx, g),
					d = data(nn_idx, g),
					thr = (thr_v.size() == G) ? thr_v[g] : thr_v[0];
				
				double
					uncertain0 = (uncertain_m.nrow() == G) ? uncertain_m(g, 0) : uncertain_m(0, 0),
					uncertain1 = (uncertain_m.nrow() == G) ? uncertain_m(g, 1) : uncertain_m(0, 1),
					missing0 = (missing_m.nrow() == G) ? missing_m(g, 0) : missing_m(0, 0),
					missing1 = (missing_m.nrow() == G) ? missing_m(g, 1) : missing_m(0, 1);
					
				x *= censor_pair(c, d, sigma, kt, thr, uncertain0, uncertain1, missing0, missing1);
				
				//Rcout << "after  i=" << row_idx << " (" << c << "), j=" << nn_idx << " (" << d << "): " << x << "\n\n";
			}
			
			triplets.push_back(T(row_idx, nn_idx, x));
		}
		if (row_idx%1000 == 0)
			callback(row_idx+1);
	}
	callback(n);
	
	Eigen::SparseMatrix<double> trans_p(n, n);
	trans_p.setFromTriplets(triplets.begin(), triplets.end());
	return trans_p;
}

// [[Rcpp::export]]
NumericMatrix predict_censoring_impl(
		NumericMatrix data,
		NumericMatrix data2,
		double thr,
		NumericVector uncertain,
		NumericVector missing,
		double sigma
) {
	int
		n = data.nrow(),
		n2 = data2.nrow(),
		G = data.ncol();
	
	if (G != data2.ncol()) stop("data and data2 must have same number of variables");
	if (uncertain.size() != 2) stop("uncertain has to be of length 2");
	if (missing.size() != 2) stop("missing has to be of length 2");
	
	//the sparse matrix has about k·n elements (or k·(n+1)?)
	NumericMatrix trans_p(n2, n);
	
	double
		kt = pow(sigma, 2),
		uncertain0 = uncertain[0] - sigma,
		uncertain1 = uncertain[1] + sigma,
		missing0 = missing[0] -= sigma,
		missing1 = missing[1] += sigma;
	
	for (int i=0; i<n; i++) {
		for (int j=0; j<n2; j++) {
			//Rcout << "i=" << row_idx << ", j=" << nn_idx << '\n';
			double x = 1.;
			for (int g=0; g<G; g++) {
				x *= censor_pair(data2(j, g), data(i, g), sigma, kt, thr, uncertain0, uncertain1, missing0, missing1);
				//Rcout << "after  i=" << row_idx << " (" << c << "), j=" << nn_idx << " (" << d << "): " << x << "\n\n";
			}
			trans_p(j, i) = x;
		}
	}
	
	return trans_p;
}
