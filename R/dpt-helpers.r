#' Find a random root cell index
#' 
#' Finds a cell that has the maximum DPT distance from a randomly selected one.
#' 
#' @param dm_or_dpt  A \code{\link{DiffusionMap}} or \code{\link{DPT}} object
#' 
#' @return A cell index
#' 
#' @examples
#' data(guo)
#' dm <- DiffusionMap(guo)
#' random_root(dm)
#' 
#' @export
random_root <- function(dm_or_dpt) {
	dpt <- dummy_dpt(dm_or_dpt)
	random_idx <- sample.int(nrow(dpt), 1L)
	which.max(dpt[random_idx, ])
}


#' Find tips in a DiffusionMap object
#' 
#' @param dm_or_dpt  A \code{\link{DiffusionMap}} or \code{\link{DPT}} object
#' @param root       Root cell index from which to find tips. (default: random)
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
find_tips <- function(dm_or_dpt, root = random_root(dm_or_dpt)) {
	dpt <- dummy_dpt(dm_or_dpt)
	tipstats(dpt, seq_len(nrow(dpt)), root)$tips
}

tipstats <- function(dpt, cells, tips) {
	x <- tips[[1L]]
	dx <- dpt[x, cells]
	y <- if (length(tips) >= 2L) tips[[2L]] else cells[[which.max(dx)]]
	dy <- dpt[y, cells]
	z <- if (length(tips) == 3L) tips[[3L]] else cells[[which.max(dx + dy)]]
	
	list(
		tips = c(x, y, z),
		dx = dx, dy = dy,
		g = max(dx + dy) / min(dx + dy))
}