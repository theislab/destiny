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
List knn_cross_impl(const NumericMatrix data, const NumericMatrix query, const size_t k, const size_t skip_self = 0) {
	if (data.ncol() != query.ncol())
		stop("data and query need the same number of features");
	const size_t nsmp_data = data.nrow();
	const size_t nsmp_query = query.nrow();
	
	std::vector<IndexedPoint<Distance>> points;
	points.reserve(nsmp_data);
	for (size_t s=0; s<nsmp_data; s++) {
		points.push_back(IndexedPoint<Distance>(data(s, _), s));
	}
	
	double max = 1e10;//std::numeric_limits<double>::infinity();
	CoverTree<IndexedPoint<Distance>> tree(max, points);
	
	IntegerMatrix index(nsmp_query, k);
	NumericMatrix dists(nsmp_query, k);
	std::fill(index.begin(), index.end(), NA_INTEGER);
	std::fill(dists.begin(), dists.end(), NA_REAL);
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(2 * nsmp_query * k);
	
	typedef std::pair<double, IndexedPoint<Distance>> DP;
	for (size_t s=0; s<nsmp_query; s++) {
		const IndexedPoint<Distance> p(query(s, _), s);
		const std::vector<DP> nns = tree.kNearestNeighborDists(p, k + skip_self);
		//Rcout << "#nn:" << nns.size() << ", k:" << k << std::endl;
		typename std::vector<DP>::const_iterator nn;
		for (nn = nns.begin() + skip_self; nn<nns.end(); nn++) {
			//Rcout << "dist:" << nn->first << ", idx:" << nn->second.idx() << ", vec:" << nn->second.vec() << std::endl;
			const size_t n = nn - nns.begin() - skip_self;
			const double dist = nn->first;
			const size_t nn_idx = nn->second.idx();
			index(s, n) = nn_idx + 1;  // R index
			dists(s, n) = dist;
			triplets.push_back(T(s, nn_idx, dist));
		}
	}
	
	Eigen::SparseMatrix<double> dist_mat(nsmp_query, nsmp_data);
	dist_mat.setFromTriplets(triplets.begin(), triplets.end());
	
	List ret;
	ret["index"] = index;
	ret["dist"]  = dists;
	ret["dist_mat"] = dist_mat;
	return ret;
}


template<class Distance>
List knn_impl(const NumericMatrix data, const size_t k) {
	return knn_cross_impl<Distance>(data, data, k, 1);
}


// [[Rcpp::export]]
List knn_cross(const NumericMatrix data, const NumericMatrix query, const size_t k, const std::string distance = "euclidean") {
	if (distance == "euclidean") {
		return knn_cross_impl<EuclideanDistance>(data, query, k);
	} else if (distance == "cosine") {
		return knn_cross_impl<CosineDistance>(data, query, k);
	} else if (distance == "rankcor") {
		NumericMatrix data_rank  = NumericMatrix(data.nrow(),  data.ncol());
		NumericMatrix query_rank = NumericMatrix(query.nrow(), query.ncol());
		for (int r=0; r<data_rank.nrow(); r++) {
			data_rank(r, _) = rank(data(r, _));
		}
		for (int r=0; r<query_rank.nrow(); r++) {
			query_rank(r, _) = rank(data(r, _));
		}
		return knn_cross_impl<CosineDistance>(data_rank, query_rank, k);
	} else stop("Unknown distance specified");
}


// [[Rcpp::export]]
List knn_asym(const NumericMatrix data, const size_t k, const std::string distance = "euclidean") {
	if (distance == "euclidean") {
		return knn_impl<EuclideanDistance>(data, k);
	} else if (distance == "cosine") {
		return knn_impl<CosineDistance>(data, k);
	} else if (distance == "rankcor") {
		NumericMatrix data_rank = NumericMatrix(data.nrow(), data.ncol());
		for (int r=0; r<data_rank.nrow(); r++) {
			data_rank(r, _) = rank(data(r, _));
		}
		return knn_impl<CosineDistance>(data_rank, k);
	} else stop("Unknown distance specified");
}
