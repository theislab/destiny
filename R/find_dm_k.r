#' Find a suitable k
#' 
#' The \code{k} parameter for the k nearest neighbors used in \link{DiffusionMap} should be as big as possible while
#' still being computationally feasible. This function approximates it depending on the size of the dataset \code{n}.
#' 
#' @param n      Number of possible neighbors (nrow(dataset) - 1)
#' @param min_k  Minimum number of neighbors. Will be chosen for \eqn{n \ge big}
#' @param small  Number of neighbors considered small. If/where \eqn{n \le small}, n itself will be returned.
#' @param big    Number of neighbors considered big. If/where \eqn{n \ge big}, \code{min_k} will be returned.
#' 
#' @return A vector of the same length as \code{n} that contains suitable \code{k} values for the respective \code{n}
#' 
#' @examples
#' curve(find_dm_k(n),     0, 13000, xname = 'n')
#' curve(find_dm_k(n) / n, 0, 13000, xname = 'n')
#' @export
find_dm_k <- function(n, min_k = 100L, small = 1000L, big = 10000L) {
	stopifnot(small < big)
	if (is.null(n)) return(NULL)
	
	k <- rep(NA_integer_, length(n))
	k[small >= n] <- n[small >= n]
	k[n >= big]   <- min_k
	
	rest <- !is.na(n) & small < n & n < big
	
	n_shifted <- (n[rest] - small) / (big - small)     # linear transf [small, big] -> [0, 1]
	k_shifted <- (cos(n_shifted * pi) + 1) / 2         # ease function [0, 1] -> [1, 0]
	k_rest    <- min_k + k_shifted * (n[rest] - min_k) # linear transf [0, 1] -> [min_k, n]
	
	k[rest] <- as.integer(round(k_rest))
	
	k
}