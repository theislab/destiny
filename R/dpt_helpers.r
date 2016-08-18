#' Find a random root cell index
#' 
#' Finds a cell index whose cell has the maximum DPT distance from a randomly selected one.
#' 
#' @param dm  A \code{\link{DiffusionMap}} object
#' 
#' @export
random_root <- function(dm) which.max(dpt_to_cell(dm, sample.int(length(dm@d_norm), 1)))


#' Find tips in a DiffusionMap object
#' 
#' @param dm    A \code{\link{DiffusionMap}} object
#' @param root  Root cell index from which to find tips. (default: random)
#' 
#' @return An integer vector of length 3
#' 
#' @export
find_tips <- function(dm, root = random_root(dm)) {
	x <- root
	dx <- dpt_to_cell(dm, x)
	y <- which.max(dx)
	dy <- dpt_to_cell(dm, y)
	z <- which.max(dx + dy)
	
	c(x, y, z)
}