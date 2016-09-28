#' Fast eigen decomposition using ARPACK
#'
#' @param M       A matrix (e.g. from the Matrix package)
#' @param n_eigs  Number of eigenvectors to return
#' @param sym     TRUE if M is symmetric
#' 
#' @return n eigenvectors of the transition matrix
#' 
#' @examples
#' eig_decomp(cbind(c(1,-1), c(-1,1)), 2)
#' 
#' @importFrom Matrix isSymmetric
#' @importFrom igraph arpack
#' @export
eig_decomp <- function(M, n_eigs, sym = isSymmetric(M)) {
	n <- nrow(M)
	f <- function(x, A = NULL) as.matrix(A %*% x)
	wh <- if (sym) 'LA' else 'LM'
	#constraints: n >= ncv > nev
	ar <- arpack(f, extra = M, sym = sym, options = list(
		which = wh, n = n, ncv = min(n, 4*n_eigs), nev = n_eigs + 1))
	if (!sym) {
		ar$vectors <- Re(ar$vectors)
		ar$values  <- Re(ar$values)
	}
	ar
}
