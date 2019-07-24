#' kNN search
#' 
#' Approximate k nearest neighbor search with flexible distance function.
#' 
#' @param data      Data matrix
#' @param query     Query matrix. Leave it out to use \code{data} as query
#' @param k         Number of nearest neighbors
#' @param ...       Parameters passed to \code{\link[RcppHNSW]{hnsw_knn}}
#' @param distance  Distance metric to use. Allowed measures: Euclidean distance (default), cosine distance (\eqn{1-corr(c_1, c_2)}) or rank correlation distance (\eqn{1-corr(rank(c_1), rank(c_2))})
#' @param sym       Return a symmetric matrix (as long as query is NULL)?
#' 
#' @return A \code{\link{list}} with the entries:
#' \describe{
#'   \item{\code{index}}{A \eqn{nrow(data) \times k} \link{integer} \link{matrix} containing the indices of the k nearest neighbors for each cell.}
#'   \item{\code{dist}}{A \eqn{nrow(data) \times k} \link{double} \link{matrix} containing the distances to the k nearest neighbors for each cell.}
#'   \item{\code{dist_mat}}{
#'     A \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} if \code{sym == TRUE},
#'     else a \code{\link[Matrix:dsCMatrix-class]{dsCMatrix}} (\eqn{nrow(query) \times nrow(data)}).
#'     Any zero in the matrix (except for the diagonal) indicates that the cells in the corresponding pair are close neighbors.
#'   }
#' }
#' 
#' @rdname knn
#' @importFrom RcppHNSW hnsw_build hnsw_knn hnsw_search
#' @export
find_knn <- function(
	data, k,
	...,
	query = NULL,
	distance = c('euclidean', 'cosine', 'rankcor', 'l2'),
	sym = TRUE,
	verbose = FALSE
) {
	p <- utils::modifyList(formals(RcppHNSW::hnsw_knn), list(...))
	if (!is.double(data)) {
		warning('find_knn does not yet support sparse matrices, converting data to a dense matrix.')
		data <- as.matrix(data)
	}
	distance <- match.arg(distance)
	
	if (distance == 'rankcor') {
		# TODO: rank_mat only works on dense matrices
		distance <- 'cosine'
		data <- rank_mat(data)
		if (!is.null(query)) query <- rank_mat(query)
	}
	
	if (is.null(query)) {
		knn <- hnsw_knn(data, k + 1L, distance, M = p$M, ef_construction = p$ef_construction, ef = p$ef, verbose = verbose)
		knn$idx  <- knn$idx[ , -1, drop = FALSE]
		knn$dist <- knn$dist[, -1, drop = FALSE]
	} else {
		index <- hnsw_build(data, distance, M = p$M, ef = p$ef_construction, verbose = verbose)
		knn <- hnsw_search(query, index, k, ef = p$ef, verbose = verbose)
	}
	names(knn)[[1L]] <- 'index'  # idx -> index
	# R matrices are column-major, so as.vector(m) == c(m[, 1], m[, 2], ...)
	knn$dist_mat <- sparseMatrix(
		rep(seq_len(nrow(knn$index)), k),
		as.vector(knn$index),
		x = as.vector(knn$dist)
	)
	if (is.null(query)) {
		if (sym) knn$dist_mat <- symmetricise(knn$dist_mat)
		nms <- rownames(data)
	} else {
		nms <- rownames(query)
	}
	rownames(knn$dist_mat) <- rownames(knn$index) <- rownames(knn$dist) <- nms
	colnames(knn$dist_mat) <- rownames(data)
	knn
}


# (double generic columnsparse to ... symmetric ...: dgCMatrix -> dsCMatrix)
# retain all differences fully. symmpart halves them in the case of trans_p[i,j] == 0 && trans_p[j,i] > 0
# TODO: could be more efficient
symmetricise <- function(dist_asym) {
	dist_sym <- symmpart(dist_asym) + abs(forceSymmetric(skewpart(dist_asym), 'U'))
	as(dist_sym, 'symmetricMatrix')
}
