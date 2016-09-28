#' Find a random root cell index
#' 
#' Finds a cell that has the maximum DPT distance from a randomly selected one.
#' 
#' @param dm  A \code{\link{DiffusionMap}} object
#' 
#' @return A cell index
#' 
#' @examples
#' data(guo)
#' dm <- DiffusionMap(guo)
#' random_root(dm)
#' 
#' @export
random_root <- function(dm) {
	acc <- accumulated_transitions(dm)
	random_idx <- sample.int(length(dm@d_norm), 1L)
	which.max(dpt_to_cell(acc, random_idx))
}


#' Find tips in a DiffusionMap object
#' 
#' @param dm    A \code{\link{DiffusionMap}} object
#' @param root  Root cell index from which to find tips. (default: random)
#' 
#' @return An integer vector of length 3
#' 
#' @examples
#' data(guo)
#' dm <- DiffusionMap(guo)
#' is_tip <- l_which(find_tips(dm), len = ncol(dataset(dm)))
#' plot(dm, col = factor(is_tip))
#' 
#' @export
find_tips <- function(dm, root = random_root(dm))
	tipstats(accumulated_transitions(dm), root)$tips

tipstats <- function(acc_trans, tips) {
	x <- tips[[1L]]
	dx <- dpt_to_cell(acc_trans, x)
	y <- if (length(tips) >= 2L) tips[[2L]] else which.max(dx)
	dy <- dpt_to_cell(acc_trans, y)
	z <- if (length(tips) == 3L) tips[[3L]] else which.max(dx + dy)
	
	list(
		tips = c(x, y, z),
		dx = dx, dy = dy,
		g = max(dx + dy) / min(dx + dy))
}