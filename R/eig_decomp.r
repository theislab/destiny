#' Fast eigen decomposition using \code{\link[RSpectra]{eigs}}
#' 
#' By default uses a random initialization vector that you can make deterministic using
#' \code{\link[base]{set.seed}} or override by specifying \code{opts = list(initvec = ...)}.
#'
#' @param M       A matrix (e.g. from the Matrix package) or
#'                a function (see \code{\link[RSpectra]{eigs}}).
#' @param n_eigs  Number of eigenvectors to return.
#' @param sym     defunct and ignored.
#' @param ...     Passed to \code{\link[RSpectra]{eigs}}.
#' 
#' @return see \code{\link[RSpectra]{eigs}}.
#' 
#' @examples
#' eig_decomp(cbind(c(1,-1), c(-1,1)), 2)
#' 
#' @importFrom Matrix isSymmetric
#' @importFrom RSpectra eigs eigs_sym
#' @export
eig_decomp <- function(M, n_eigs, sym, ..., opts = list()) {
	if (!('initvec' %in% names(opts)))
		opts$initvec <- runif(nrow(M)) - .5
	# eigs cannot handle symmetricMatrix & sparseMatrix yet
	# TODO: low-effort. We can copy the memory and use the `lower = T/F` arg instead
	if (is(M, 'dsCMatrix')) M <- as(M, 'dgCMatrix')
	if (is(M, 'dsRMatrix')) M <- as(M, 'dgRMatrix')
	eigs(M, n_eigs, ..., opts = opts)
}
