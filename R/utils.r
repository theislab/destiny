stopifsmall <- function(max_dist) {
	if (max_dist < .Machine$double.eps)
		stop(sprintf(
			'The supplied sigma is not large enough. Please select a larger one.
			find_sigmas(data) should return one with the right order of magnitude. (max dist. is %.3e)',
			max_dist))
}


#' @importFrom utils flush.console
verbose_timing <- function(verbose, msg, expr) {
	if (verbose) {
		cat(sprintf('%s...', msg))
		flush.console()
		dif <- system.time({
			r <- force(expr)
		})
		cat(sprintf('...done. Time: %.2fs\n', dif[['elapsed']]))
		flush.console()
		r
	} else expr
}


#' @importFrom Matrix Diagonal
#' @importMethodsFrom Matrix solve
accumulated_transitions <- function(dm) {
	if (!is.null(dm@data_env$propagations)) {  # compat
		dm@data_env$accumulated_transitions <- dm@data_env$propagations
		rm('propagations', envir = dm@data_env)
	}
	
	if (is.null(dm@data_env$accumulated_transitions)) {
		if (is.null(dm@transitions))
			stop('DiffusionMap was created with suppress_dpt = TRUE')
		
		n <- length(dm@d_norm)
		
		phi0 <- dm@d_norm / sqrt(sum(dm@d_norm ^ 2))
		inv <- solve(Diagonal(n) - dm@transitions + phi0 %*% t(phi0))
		dm@data_env$accumulated_transitions <- inv - Diagonal(n)
	}
	
	dm@data_env$accumulated_transitions
}


hasattr <- function(x, which) !is.null(attr(x, which, exact = TRUE))


flipped_dcs <- function(d, dcs) {
	if (is(d, 'DiffusionMap')) d <- eigenvectors(d)
	
	evs <- as.matrix(d[, abs(dcs)])
	evs[, dcs < 0] <- -evs[, dcs < 0]
	evs
}


rescale_mat <- function(mat, rescale) {
	if (is.list(rescale)) {
		stopifnot(setequal(dimnames(rescale), c('from', 'to')))
		rv <- apply(mat, 2L, scales::rescale, rescale$to, rescale$from)
	} else if (is.array(rescale)) {
		stopifnot(length(dim(rescale)) == 3L)
		stopifnot(ncol(mat) == ncol(rescale))
		stopifnot(dim(rescale)[[1L]] == 2L)
		stopifnot(dim(rescale)[[3L]] == 2L)
		
		col_type <- get(typeof(mat))
		rv <- vapply(seq_len(ncol(mat)), function(d) {
			scales::rescale(mat[, d], rescale['to', d, ], rescale['from', d, ])
		}, col_type(nrow(mat)))
	}
	stopifnot(all(dim(rv) == dim(mat)))
	dimnames(rv) <- dimnames(mat)
	rv
}


runs <- function(vec) {
	enc <- rle(vec)
	enc$values <- make.unique(enc$values, '_')
	inverse.rle(enc)
}


upper.tri.sparse <- function(x, diag = FALSE) {
	# Works just like upper.tri() but doesn't forcibly coerce large 'sparseMatrix' back to 'matrix'
	if (diag)
		row(x) <= col(x)
	else row(x) < col(x)
}


#' @importFrom igraph graph_from_adjacency_matrix membership cluster_louvain
get_louvain_clusters <- function(transitions) {
	graph <- graph_from_adjacency_matrix(transitions, 'undirected', weighted = TRUE)
	as.integer(unclass(membership(cluster_louvain(graph))))
}
