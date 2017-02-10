#' @include diffusionmap.r
NULL

#' Gene relevances for entire data set
#' 
#' The relevance map is cached insided of the \code{\link{DiffusionMap}}.
#' 
#' @param coords  A \code{\link{DiffusionMap}} object or a cells \eqn{\times} dims \code{\link{matrix}}.
#' @param exprs   An cells \eqn{\times} genes \code{\link{matrix}}. Only provide if \code{coords} is no \code{\link{DiffusionMap}}.
#' @param ...     Ignored.
#' @param k       Number of nearest neighbors to use
#' @param dims    Index into columns of \code{coord}
#' 
#' @return A \code{GeneRelevance} object:
#' 
#' @slot coords         A cells \eqn{\times} dims \code{\link{matrix}} of coordinates (e.g. diffusion components), reduced to the dimensions passed as \code{dims}.
#' @slot exprs          A cells \eqn{\times} genes matrix of expressions.
#' @slot partials       Array of partial derivatives wrt to considered dimensions in reduced space.
#'                      (genes \eqn{\times} cells \eqn{\times} dimensions)
#' @slot partials_norm  Matrix with norm of aforementioned derivatives. (n\_genes \eqn{\times} cells)
#' @slot nn_index       Matrix of k nearest neighbor indices. (cells \eqn{\times} k)
#' 
#' @seealso \link{Gene Relevance plotting}: \code{plot_gradient_map}/\code{plot_gene_relevance}
#' 
#' @aliases gene_relevance
#' @name Gene Relevance
#' @export
setClass('GeneRelevance', slots = c(
	coords = 'matrix',
	exprs = 'matrix',
	partials = 'array',
	partials_norm = 'matrix',
	nn_index = 'matrix'))

#' @name Gene Relevance
#' @export
setGeneric('gene_relevance', function(coords, exprs, ..., k = 20, dims = 1:2, verbose = FALSE) standardGeneric('gene_relevance'))

#' @name Gene Relevance
#' @export
setMethod('gene_relevance', c('DiffusionMap', 'missing'), function(coords, exprs, ..., k = 20, dims = 1:2, verbose = FALSE) {
	dm <- coords
	if (!is.null(dm@data_env$relevance_map))
		return(dm@data_env$relevance_map)
	coords <- eigenvectors(dm)
	exprs <- extract_doublematrix(dataset(dm))
	dm@data_env$relevance_map <- gene_relevance_impl(coords, exprs, k = k, dims = dims, verbose = verbose)
	dm@data_env$relevance_map
})

#' @name Gene Relevance
#' @export
setMethod('gene_relevance', c('matrix', 'matrix'), function(coords, exprs, ..., k = 20, dims = 1:2, verbose = FALSE) {
	gene_relevance_impl(coords, exprs, k = k, dims = dims, verbose = verbose)
})

#' @importFrom FNN get.knn
#' @importFrom Biobase rowMedians
gene_relevance_impl <- function(coords, exprs, ..., k, dims, verbose) {
	if (is.null(colnames(exprs))) stop('The expression matrix columns need to be named but are NULL')
	coords <- coords[, dims]
	nn_index <- get.knn(exprs, k, algorithm = 'cover_tree')$nn.index
	
	k <- ncol(nn_index)
	n_genes <- ncol(exprs)
	n_cells <- nrow(coords)
	n_dims <- ncol(coords) # Diffusion components to compute partial derivatives for
	partials <- array(NA, dim = c(n_genes, n_cells, n_dims))
	
	if (verbose) cat('Calculating expression gradient\n')
	gradient_exprs <- apply(exprs, 2L, function(expr_gene) {
		# Compute change in expression
		# Do not compute if reference is zero, could be drop-out
		expr_masked <- expr_gene
		expr_masked[expr_masked == 0] <- NA
		gradient_expr <- apply(nn_index, 2, function(nn) expr_gene[nn] - expr_masked)
		gradient_expr[gradient_expr == 0] <- NA  # Cannot evaluate partial
		stopifnot(length(dim(gradient_expr)) == 2L)
		gradient_expr
	})
	# apply only handles returning vectors, so we have to reshape the return value
	dim(gradient_exprs) <- c(n_cells, k, n_genes)
	
	for (d in seq_len(n_dims)) {
		# Compute partial derivatives in direction of current dimension
		
		if (verbose) cat('Calculating partial derivatives of dimension ', d, '/', n_dims, '\n')
		# We could optionaly add normalization by max(coords[, d]) - min(coords[, d])
		gradient_coord <- apply(nn_index, 2L, function(nn) coords[nn, d] - coords[, d])
		
		partials[, , d] <- apply(gradient_exprs, 3L, function(grad_gene_exprs) {
			# Compute median of difference quotients to NN
			difference_quotients <- gradient_coord / grad_gene_exprs
			# Only compute gradient if at least two observations are present!
			stable_cells <- rowSums(!is.na(difference_quotients)) >= 2L
			ifelse(stable_cells, rowMedians(difference_quotients), NA)
		})
	}
	
	# Compute norm over partial derivates: Frobenius
	partials_norm <- apply(partials, c(1, 2), function(z) sqrt(sum(z^2, na.rm = TRUE)))
	
	# Find outlier cells: Not in NN of more than 1 other cell
	# Remove these as they tend to receive very large norms
	#outliers <- sapply(seq_len(n_cells), function(cell) sum(nn_index == cell) > 1)
	#partials_norm[, outliers] <- NA
	
	# Prepare output
	rownames(partials_norm) <- colnames(exprs)
	dimnames(partials)[[1]] <- colnames(exprs)
	dimnames(partials)[[3]] <- if (is.character(dims)) dims else colnames(coords)
	new('GeneRelevance',
		coords = coords,
		exprs = exprs,
		partials = partials,
		partials_norm = partials_norm,
		nn_index = nn_index)
}
