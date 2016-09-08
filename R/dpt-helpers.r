#' Find a random root cell index
#' 
#' Finds a cell index whose cell has the maximum DPT distance from a randomly selected one.
#' 
#' @param dm  A \code{\link{DiffusionMap}} object
#' 
#' @export
random_root <- function(dm) {
	propagations <- propagation_matrix(dm)
	random_idx <- sample.int(length(dm@d_norm), 1L)
	which.max(dpt_to_cell(propagations, random_idx))
}


#' Find tips in a DiffusionMap object
#' 
#' @param dm    A \code{\link{DiffusionMap}} object
#' @param root  Root cell index from which to find tips. (default: random)
#' 
#' @return An integer vector of length 3
#' 
#' @export
find_tips <- function(dm, root = random_root(dm))
	tipstats(propagation_matrix(dm), root)$tips

tipstats <- function(propagations, tips) {
	x <- tips[[1L]]
	dx <- dpt_to_cell(propagations, x)
	y <- if (length(tips) >= 2L) tips[[2L]] else which.max(dx)
	dy <- dpt_to_cell(propagations, y)
	z <- if (length(tips) == 3L) tips[[3L]] else which.max(dx + dy)
	
	list(
		tips = c(x, y, z),
		dx = dx, dy = dy,
		g = max(dx + dy) / min(dx + dy))
}