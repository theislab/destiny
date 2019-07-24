#' @include diffusionmap.r s4-unions.r
NULL

#' Gene relevances for entire data set
#' 
#' The relevance map is cached insided of the \code{\link{DiffusionMap}}.
#' 
#' @param coords    A \code{\link{DiffusionMap}} object or a cells \eqn{\times} dims \code{\link{matrix}}.
#' @param exprs     An cells \eqn{\times} genes \code{\link{matrix}}. Only provide if \code{coords} is no \code{\link{DiffusionMap}}.
#' @param ...       If no \code{\link{DiffusionMap}} is provided, a vector of \code{weights} (of the same length as \code{dims}) can be provided.
#' @param k         Number of nearest neighbors to use
#' @param dims      Index into columns of \code{coord}
#' @param distance  Distance measure to use for the nearest neighbor search.
#' @param smooth    Smoothing parameters \code{c(window, alpha)} (see \code{\link[smoother]{smth.gaussian}}).
#'                  Alternatively \code{\link{TRUE}} to use the \link[smoother]{smoother} \link[smoother:smth.options]{defaults}
#'                  or \code{\link{FALSE}} to skip smoothing,
#' @param verbose   If TRUE, log additional info to the console
#' 
#' @return A \code{GeneRelevance} object:
#' 
#' @slot coords         A cells \eqn{\times} dims \code{\link{matrix}} or \code{\link[Matrix:sparseMatrix-class]{sparseMatrix}}
#'                      of coordinates (e.g. diffusion components), reduced to the dimensions passed as \code{dims}
#' @slot exprs          A cells \eqn{\times} genes matrix of expressions
#' @slot partials       Array of partial derivatives wrt to considered dimensions in reduced space
#'                      (genes \eqn{\times} cells \eqn{\times} dimensions)
#' @slot partials_norm  Matrix with norm of aforementioned derivatives. (n\_genes \eqn{\times} cells)
#' @slot nn_index       Matrix of k nearest neighbor indices. (cells \eqn{\times} k)
#' @slot dims           Column index for plotted dimensions. Can \code{\link{character}}, \code{\link{numeric}} or \code{\link{logical}}
#' @slot distance       Distance measure used in the nearest neighbor search. See \code{\link{find_knn}}
#' @slot smooth_window  Smoothing window used (see \code{\link[smoother]{smth.gaussian}})
#' @slot smooth_alpha   Smoothing kernel width used (see \code{\link[smoother]{smth.gaussian}})
#' 
#' @seealso \link{Gene Relevance methods}, \link{Gene Relevance plotting}: \code{plot_differential_map}/\code{plot_gene_relevance}
#' 
#' @examples
#' data(guo_norm)
#' dm <- DiffusionMap(guo_norm)
#' gr <- gene_relevance(dm)
#' 
#' m <- t(Biobase::exprs(guo_norm))
#' gr_pca <- gene_relevance(prcomp(m)$x, m)
#' # now plot them!
#' 
#' @rdname Gene-Relevance
#' @export
setClass('GeneRelevance', slots = c(
	coords = 'matrix',
	exprs = 'dMatrixOrMatrix',
	partials = 'array',
	partials_norm = 'matrix',
	nn_index = 'matrix',  # k = ncol(nn_index)
	dims = 'ColIndex',
	distance = 'character',
	smooth_window = 'numeric',
	smooth_alpha = 'numeric'))

#' @rdname Gene-Relevance
#' @export
setGeneric('gene_relevance', function(coords, exprs, ..., k = 20L, dims = 1:2, distance = NULL, smooth = TRUE, remove_outliers = FALSE, verbose = FALSE) standardGeneric('gene_relevance'))

#' @importFrom Biobase updateObject
#' @rdname Gene-Relevance
#' @export
setMethod('gene_relevance', c('DiffusionMap', 'missing'), function(coords, exprs, ..., k, dims, distance, smooth, remove_outliers, verbose) {
	dm <- coords
	relevance_map <- updateObject(dm@data_env$relevance_map)
	smooth <- get_smoothing(smooth)
	if (is.null(relevance_map) ||
		ncol(relevance_map@nn_index) != k ||
		!identical(relevance_map@dims, dims) ||
		!identical(relevance_map@smooth_window, smooth[[1L]]) ||
		!identical(relevance_map@smooth_alpha,  smooth[[2L]])
	) {
		coords <- eigenvectors(dm)
		exprs <- dataset_extract_doublematrix(dataset(dm))
		weights <- eigenvalues(dm)[dims]
		if (is.null(distance)) distance <- dm@distance
		else if (!identical(distance, dm@distance)) stop('the specified distance ', distance,' is not the same as the one used for the diffusion map: ', dm@distance)
		relevance_map <- gene_relevance_impl(coords, exprs, ..., k = k, dims = dims, distance = distance, smooth = smooth, remove_outliers = remove_outliers, verbose = verbose, weights = weights)
		dm@data_env$relevance_map <- relevance_map
	} else stopifparams(...)
	relevance_map
})

#' @rdname Gene-Relevance
#' @export
setMethod('gene_relevance', c('matrix', 'matrix'), function(coords, exprs, ..., k, dims, distance, smooth, remove_outliers, verbose) {
	gene_relevance_impl(coords, exprs, ..., k = k, distance = distance, smooth = smooth, dims = dims, remove_outliers = remove_outliers, verbose = verbose)
})

#' @importFrom Biobase rowMedians
gene_relevance_impl <- function(coords, exprs, ..., k, dims, distance, smooth, remove_outliers, verbose, weights = 1) {
	stopifparams(...)
	distance <- match.arg(distance, c('euclidean', 'cosine', 'rankcor', 'l2'))
	coords_used <- coords[, dims, drop = FALSE]
	n_dims <- ncol(coords_used)
	if (length(weights) == 1L) weights <- rep(weights, n_dims)
	smooth <- get_smoothing(smooth)
	
	if (is.null(colnames(exprs))) stop('The expression matrix columns need to be named but are NULL')
	if (n_dims != length(weights)) stop(n_dims, ' dimensions, but ', length(weights), ' weights were provided')
	
	nn_index <- find_knn(exprs, k, distance = distance)$index
	
	k <- ncol(nn_index)
	n_cells <- nrow(coords_used)
	n_genes <- ncol(exprs)
	partials <- array(
		NA,
		dim = c(n_cells, n_genes, n_dims),
		dimnames = list(rownames(exprs), colnames(exprs), if (is.character(dims)) dims else colnames(coords_used)))
	
	# a very small value to subtract from the differential
	small <- min(exprs[exprs != 0]) / length(exprs[exprs == 0])
	if (verbose) cat('Calculating expression differential\n')
	gene_differential <- function(expr_gene) {
		# Compute change in expression
		# Do not compute if reference is zero, could be drop-out
		expr_masked <- expr_gene
		expr_masked[expr_masked == 0] <- small
		differential_expr <- apply(nn_index, 2, function(nn) expr_gene[nn] - expr_masked)
		differential_expr[differential_expr == 0] <- NA  # Cannot evaluate partial
		#stopifnot(identical(dim(differential_expr), c(n_cells, k)))
		differential_expr
	}
	differential_exprs <- apply(exprs, 2L, gene_differential)
	#stopifnot(identical(dim(differential_exprs), c(n_cells * k, n_genes)))
	# apply only handles returning vectors, so we have to reshape the return value
	dim(differential_exprs) <- c(n_cells, k, n_genes)
	dimnames(differential_exprs)[[3L]] <- if (length(colnames(exprs)) > 1L) colnames(exprs) else list(colnames(exprs))
	#stopifnot(identical(gene_differential(exprs[, 1L]), differential_exprs[, , 1L]))
	
	for (d in seq_len(n_dims)) {
		# Compute partial derivatives in direction of current dimension
		
		if (verbose) cat('Calculating partial derivatives of dimension ', d, '/', n_dims, '\n')
		# We could optionaly add normalization by max(coords_used[, d]) - min(coords_used[, d])
		differential_coord <- apply(nn_index, 2L, function(nn) coords_used[nn, d] - coords_used[, d])
		
		partials_unweighted <- apply(differential_exprs, 3L, function(grad_gene_exprs) {
			# Compute median of difference quotients to NN
			difference_quotients <- grad_gene_exprs / differential_coord
			# Only compute differential if at least two observations are present!
			stable_cells <- rowSums(!is.na(difference_quotients)) >= 2L
			ifelse(stable_cells, rowMedians(difference_quotients, na.rm = TRUE), NA)
		})
		colnames(partials_unweighted) <- colnames(exprs)
		
		if (!any(is.na(smooth))) {
			order_coor <- order(coords_used[, d])
			order_orig <- order(order_coor)
			# Smooth the partials for each gene across all cells
			partials_unweighted <- apply(partials_unweighted, 2L, function(partials_gene) {
				ordered <- partials_gene[order_coor]
				smoothed <- smth.gaussian(ordered, smooth[[1L]], smooth[[2L]], tails = TRUE)
				smoothed[order_orig]
			})
			colnames(partials_unweighted) <- colnames(exprs)
		}
		
		partials[, , d] <- weights[[d]] * partials_unweighted
	}
	
	# Compute norm over partial derivates: Frobenius
	partials_norm <- apply(partials, c(1, 2), function(z) sqrt(sum(z^2, na.rm = TRUE)))
	colnames(partials_norm) <- colnames(partials)
	
	# Find outlier cells: Not in NN of more than 1 other cell
	# Remove these as they tend to receive very large norms
	if (remove_outliers) {
		outliers <- sapply(seq_len(n_cells), function(cell) sum(nn_index == cell) > 1)
		partials_norm[, outliers] <- NA
	}
	
	# Prepare output
	rownames(partials_norm) <- rownames(partials)
	colnames(partials_norm) <- colnames(partials)
	new('GeneRelevance',
		coords = coords,
		exprs = exprs,
		partials = partials,
		partials_norm = partials_norm,
		nn_index = nn_index,
		dims = dims,
		smooth_window = smooth[[1L]],
		smooth_alpha  = smooth[[2L]],
		distance = distance)
}

get_smoothing <- function(smooth) {
	if (isTRUE(smooth)) c(getOption('smoother.window'), getOption('smoother.gaussianwindow.alpha'))
	else if (identical(smooth, FALSE)) c(NA_real_, NA_real_)
	else if (!is.numeric(smooth) || length(smooth) != 2L)
		stop('`smooth` needs to be TRUE, FALSE or a numeric c(window, alpha), not', capture.output(str(smooth)))
	else smooth
}
