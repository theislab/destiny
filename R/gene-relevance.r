#' @include diffusionmap.r
NULL

#' Gene relevances for entire data set
#' 
#' The relevance map is cached insided of the \code{\link{DiffusionMap}}.
#' 
#' @param dm      A \code{\link{DiffusionMap}} object.
#' @param dims    Index into columns of \code{coord}
#' @param n_proc  Number of threads for parallelisation.
#' 
#' @return A \code{GeneRelevance} object:
#' 
#' @slot partial        Array of partial derivatives wrt to considered dimensions in reduced space.
#'                      (genes \eqn{\times} cells \eqn{\times} dimensions)
#' @slot partials_norm  Matrix with norm of aforementioned derivatives. (n\_genes \eqn{\times} cells)
#' @slot nn_index       Matrix of k nearest neighbor indices. (cells \eqn{\times} k)
#' 
#' @seealso \link{Gene Relevance plotting}: \code{plot_gradient_map}/\code{plot_gene_relevance}
#' 
#' @importFrom BiocParallel register MulticoreParam SerialParam bplapply
#' @importFrom FNN get.knn
#' 
#' @aliases gene_relevance
#' @name Gene Relevance
#' @export
setClass('GeneRelevance', slots = c(
	dm = 'DiffusionMap',
	partials = 'array',
	partials_norm = 'matrix',
	nn_index = 'matrix'))

#' @name Gene Relevance
#' @export
gene_relevance <- function(dm, k = 20, dims = 1:2, n_proc = 1) {
	if (!is(dm, 'DiffusionMap')) stop('dm needs to be a DiffusionMap object')
	if (!is.null(dm@data_env$relevance_map)) return(dm@data_env$relevance_map)
	
	coords <- dm@eigenvectors
	exprs <- extract_doublematrix(dataset(dm))
	
	# Set parallelisation
	bpparam <- if (n_proc > 1) MulticoreParam(workers = n_proc) else SerialParam()
	#partials[, , dc] <- do.call(rbind, bplapply(seq_len(n_genes), BPPARAM = bpparam, FUN = function(gene) {
	
	coords <- coords[, dims]
	nn_index <- get.knn(exprs, k, algorithm = 'cover_tree')$nn.index
	
	k <- ncol(nn_index)
	n_genes <- ncol(exprs)
	n_cells <- nrow(coords)
	n_dims <- ncol(coords) # Diffusion components to compute partial derivatives for
	partials <- array(NA, dim = c(n_genes, n_cells, n_dims))
	
	gradient_exprs <- apply(exprs, 2L, function(expr_gene) {
		# Compute change in expression
		# Do not compute if reference is zero, could be drop-out
		expr_masked <- expr_gene
		expr_masked[expr_masked == 0] <- NA
		gradient_expr <- apply(nn_index, 2, function(nn) expr_gene[nn] - expr_masked)
		gradient_expr[gradient_expr == 0] <- NA  # Cannot evaluate partial
		stopifnot(length(dim(gradient_expr)) == 2)
		gradient_expr
	})
	dim(gradient_exprs) <- c(n_cells, k, n_genes)
	
	for (dc in seq_len(n_dims)) {
		# Compute partial derivatives in direction of current 
		# dimension (e.g. current DC)
		
		# We could optionaly add normalization by  max(coords[,dc]) - min(coords[,dc])
		gradient_coord <- apply(nn_index, 2, function(nn) coords[nn,dc] - coords[,dc])
		
		# Parallise over genes:
		partials[, , dc] <- apply(seq_len(n_genes), function(gene) {
			# Compute median of difference quotients to NN
			apply(gradient_coord / gradient_exprs[, , gene], 1, function(cell) {
				# Only compute gradient if at least two observations are present!
				# Otherwise very unstable
				if (sum(!is.na(cell)) > 1) median(cell, na.rm = TRUE) else NA
			})
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
	dimnames(partials)[[3]] <- if (is.character(dims)) dims else colnames(coords)[dims]
	(dm@data_env$relevance_map <- new('GeneRelevance',
		dm = dm,
		partials = partials,
		partials_norm = partials_norm,
		nn_index = nn_index))
}
