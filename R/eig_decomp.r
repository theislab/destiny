#' Fast eigen decomposition using \code{\link[RSpectra]{eigs}}
#'
#' By default uses a random initialization vector that you can make deterministic using
#' \code{\link[base]{set.seed}} or override by specifying \code{opts = list(initvec = ...)}.
#'
#' @param m       A matrix (e.g. from the Matrix package) or
#'                a function (see \code{\link[RSpectra]{eigs}}).
#' @param n_eigs  Number of eigenvectors to return.
#' @param sym     defunct and ignored.
#' @param ...     Passed to \code{\link[RSpectra]{eigs}}.
#' @param opts    Passed to \code{\link[RSpectra]{eigs}}.
#'
#' @return see \code{\link[RSpectra]{eigs}}.
#'
#' @examples
#' eig_decomp(cbind(c(1,0,-1), c(0,1,0), c(-1,0,1)), 2)
#'
#' @importFrom Matrix isSymmetric
#' @importFrom RSpectra eigs eigs_sym
#' @importFrom stats runif
#' @export
eig_decomp <- function(m, n_eigs, sym, ..., opts = list()) {
	if (!('initvec' %in% names(opts)))
		opts$initvec <- runif(nrow(m)) - .5
	# eigs cannot handle symmetricMatrix & sparseMatrix yet
	# TODO: low-effort. We can copy the memory and use the `lower = T/F` arg instead
	if (is(m, 'dsCMatrix')) m <- as(m, 'dgCMatrix')
	if (is(m, 'dsRMatrix')) m <- as(m, 'dgRMatrix')
	eigs(m, n_eigs, ..., opts = opts)
}
