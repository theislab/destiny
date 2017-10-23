#' kNN search
#' 
#' k nearest neighbor search with custom distance function.
#' 
#' @param data      Data matrix
#' @param query     Query matrix. In \code{knn} and \code{knn_asym}, query and data are identical.
#' @param k         Number of nearest neighbors.
#' @param ...       Ignored.
#' @param distance  Distance metric to use. Allowed measures: Euclidean distance (default), cosine distance (\eqn{1-corr(c_1, c_2)}) or rank correlation distance (\eqn{1-corr(rank(c_1), rank(c_2))}).
#' @param sym       Return a symmetric matrix (as long as query is NULL)?
#' 
#' @return A \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} if \code{sym == TRUE}, else a \code{\link[Matrix:dsCMatrix-class]{dsCMatrix}} (\eqn{nrow(query) \times nrow(data)}).
#' 
#' @name knn
#' @export
find_knn <- function(data, k, ..., query = NULL, distance = c('euclidean', 'cosine', 'rankcor'), sym = TRUE) {
	distance <- match.arg(distance)
	if (is.null(query)) {
		knn <- knn_asym(data, k, distance)
		if (sym) knn$dist_mat <- symmetricise(knn$dist_mat)
		knn
	} else knn_cross(data, query, k, distance)
}


# (double generic columnsparse to ... symmetric ...: dgCMatrix -> dsCMatrix)
# retain all differences fully. symmpart halves them in the case of trans_p[i,j] == 0 && trans_p[j,i] > 0
# TODO: could be more efficient
symmetricise <- function(dist_asym) {
	dist_sym <- symmpart(dist_asym) + abs(forceSymmetric(skewpart(dist_asym), 'U'))
	as(dist_sym, 'symmetricMatrix')
}
