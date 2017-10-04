#include <Rcpp.h>
#include <RcppEigen.h>
#include "Cover_Tree.h"
#include "utils.h"

using namespace Rcpp;


template<class Distance>
class IndexedPoint {
private:
	NumericVector _vec;
	size_t _idx;
public:
	IndexedPoint(NumericVector v, size_t i) : _vec(v), _idx(i) {}
	const NumericVector& vec() const { return this->_vec; }
	size_t               idx() const { return this->_idx; }
	bool operator==(const IndexedPoint<Distance>& p) {
		return is_true(all(this->_vec == p.vec()));
	};
	double distance(const IndexedPoint<Distance>& p) const {
		return Distance::distance(*this, p);
	};
};


template<class Distance>
std::ostream &operator<<(std::ostream& os, IndexedPoint<Distance> const& m) { 
	return os << "IndexedPoint(" << m.idx() << ", " << m.vec() << ")";
}


class EuclideanDistance {
public:
	static double distance(const IndexedPoint<EuclideanDistance>& p1, const IndexedPoint<EuclideanDistance>& p2) {
		const NumericVector& vec1 = p1.vec();
		const NumericVector& vec2 = p2.vec();
		
		double dist = 0;
		const size_t lim = vec1.size();
		for (size_t i = 0; i<lim; i++) {
			double d = vec1[i] - vec2[i];
			dist += d*d;
		}
		dist = sqrt(dist);
		
		return dist;
	}
};


class CosineDistance {
public:
	static double distance(const IndexedPoint<CosineDistance>& p1, const IndexedPoint<CosineDistance>& p2) {
		return 1 - cor(p1.vec(), p2.vec());
	}
};


template<class Distance>
List knn_impl(NumericMatrix data, size_t k) {
	const size_t n_samples = data.nrow();
	const size_t n_features = data.ncol();
	
	std::vector<IndexedPoint<Distance>> points;
	points.reserve(n_samples);
	for (size_t s=0; s<n_samples; s++) {
		points.push_back(IndexedPoint<Distance>(data(s, _), s));
	}
	
	double max = 1e10;//std::numeric_limits<double>::infinity();
	CoverTree<IndexedPoint<Distance>> tree(max, points);
	
	IntegerMatrix index(n_samples, k);
	NumericMatrix dists(n_samples, k);
	std::fill(index.begin(), index.end(), NA_INTEGER);
	std::fill(dists.begin(), dists.end(), NA_REAL);
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(2 * n_samples * k);
	
	typedef std::pair<double, IndexedPoint<Distance>> DP;
	for (size_t s=0; s<n_samples; s++) {
		// skip first NN, as it’s the point itself.
		const std::vector<DP> nns = tree.kNearestNeighborDists(points[s], k + 1);
		//Rcout << "#nn:" << nns.size() << ", k:" << k << std::endl;
		typename std::vector<DP>::const_iterator nn;
		for (nn = nns.begin() + 1; nn<nns.end(); nn++) {
			//Rcout << "dist:" << nn->first << ", idx:" << nn->second.idx() << ", vec:" << nn->second.vec() << std::endl;
			const size_t n = nn - nns.begin() - 1;
			const double dist = nn->first;
			index(s, n) = nn->second.idx() + 1;  // R index
			dists(s, n) = dist;
			triplets.push_back(T(s, n, dist));
			triplets.push_back(T(n, s, dist));
		}
	}
	
	Eigen::SparseMatrix<double> dist_mat(n_samples, n_samples);
	dist_mat.setFromTriplets(triplets.begin(), triplets.end());
	
	List ret;
	ret["index"] = index;
	ret["dist"]  = dists;
	ret["dist_mat"] = dist_mat;
	return ret;
}

// [[Rcpp::export]]
List knn(NumericMatrix imputed_data, size_t k, std::string distance = "euclidean") {
	if (distance == "euclidean") {
		return knn_impl<EuclideanDistance>(imputed_data, k);
	} else if (distance == "cosine") {
		return knn_impl<CosineDistance>(imputed_data, k);
	} else if (distance == "rank") {
		NumericMatrix data = NumericMatrix(imputed_data.nrow(), imputed_data.ncol());
		for (int r=0; r<data.nrow(); r++) {
			data(r, _) = rank(imputed_data(r, _));
		}
		return knn_impl<CosineDistance>(data, k);
	} else stop("Unknown distance specified");
}
