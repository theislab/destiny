#' Split a branch into three
#' 
#' Assign cells to a branch b by maximizing finite Kendall correlation of dpts of the
#' other two branches on b + finite Kendall anticorrelation of dpts on the other branches.
#' 
#' Both dpt and indb contain an entry for each cell, as this function performs the cut.
#'
#' @param dpt  Pseudo time distances for each of the three branches     (\eqn{n \times 3} matrix of double)
#' @param bid  Cell indices in ascending pseudotime distance per branch (\eqn{n \times 3} matrix of integer)
#' @param b    Index of branch to cut (1, 2, or 3)
#' 
#' @importFrom smoother smth.gaussian
#' @export
branchcut <- function(dpt, bid, b) {
	n <- nrow(bid)
	all_branches <- seq_len(3L)
	
	# sanity checks
	stopifnot(b %in% all_branches)
	stopifnot(ncol(dpt) == 3L, ncol(bid) == 3L)
	stopifnot(nrow(dpt) == n)
	stopifnot(is.double(dpt), is.integer(bid))
	
	# find cell indexes per branch 
	other <- all_branches[all_branches != b]
	b1 <- other[[1L]]
	b2 <- other[[2L]]
	
	# DPT for other branches, sorted by b3
	b3_idxs <- bid[, b]
	dpt1 <- dpt[b3_idxs, b1]
	dpt2 <- dpt[b3_idxs, b2]
	
	kcor <- vapply(seq_len(n - 1L), function(s1) {
		s2 <- s1 + 1L
		l <- seq_len(s1)
		r <- seq(s2, n)
		
		k_l <- kendall_finite_cor(dpt1[l], dpt2[l], dpt1[[s2]], dpt2[[s2]])
		k_r <- kendall_finite_cor(dpt1[r], dpt2[r], dpt1[[s1]], dpt2[[s1]])
		
		k_l/s1 - k_r/(n - s1)
	}, double(1))
	
	kcor <- smth.gaussian(kcor, 5L)
	cut <- which.max(kcor)
	
	b3_idxs[seq_len(cut)]
}


# This function organizes the cell arrays branch[[i]]
# to build Branch labels (length <- number of cells) indicating the branch each cell belongs to.
# Cells which are assigned to more than one branch in dpt as well
# as cells which are not assigned to any branch are defined as undeceided (label NA)
#' @importFrom utils combn
organize_branches <- function(branch) {
	n <- do.call(max, branch)
	
	intersect_branches <- function(bs) intersect(branch[[bs[[1]]]], branch[[bs[[2]]]])
	branch_intersections <- lapply(combn(3L, 2L, simplify = FALSE), intersect_branches)
	inters <- Reduce(union, branch_intersections, integer())
	
	branch <- lapply(branch, function(b) setdiff(b, inters))
	branch_nums <- seq_along(branch)  # TODO: change
	
	branches <- rep(NA_integer_, n)
	for (b in seq_along(branch)) {
		branches[branch[[b]]] <- branch_nums[[b]]
	}
	branches
}


# given two orderings b1 and b2, compute the delta (e.i. finite) Kendall
# correlation of adding a new cell with bnew1, bnew2 to the orderings.
kendall_finite_cor <- function(b1, b2, new1, new2) {
	b11 <- numeric(length(b1))
	b11[b1 >= new1] <- 1
	b11[b1 <  new1] <- -1
	
	b22 <- numeric(length(b2))
	b22[b2 >= new2] <- 1
	b22[b2 <  new2] <- -1
	
	b11 %*% b22
}
