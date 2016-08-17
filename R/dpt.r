#' Diffusion Pseudo Time
#'
#' Create pseudotime ordering and assigns cell to one of three branches
#' 
#' @param dm           A \code{\link{DiffusionMap}} object
#' @param branching    Detect a branching? (\code{TRUE} or \code{FALSE})
#' @param tips         Tip cell indices for each branch (integer vector of length 1 or 3)
#' @param root         The root index from which to calculate the DPTs (integer of length 1)
#' 
#' @details
#' All parameters are optional, but at least \code{branching} or \code{tips} has to be there.
#' If unspecified: \describe{
#'  \item{\code{root}}{will be the furthest cell from a random one}
#' 	\item{\code{branching}}{will be \code{TRUE} if multiple \code{tips} are specified}
#'  \item{\code{tips}}{will be \code{root} and (\code{if (branching)}) the most distant cells from it}
#' }
#' 
#' @return A \code{\link[base]{data.frame}} with the rows: \describe{
#' 	\item{\code{Branch}}{Branch labels for each cell, 1,2,3 or NA for undeceided}
#'  \item{\code{DPT}}{Diffusion pseudotime in respect to the root cell}
#'  \item{(\code{if (branching)}) \code{DPT.1}, \code{DPT.2}}{Diffusion pseudotime in respect to the other tips}
#' }
#' 
#' @export
dpt <- function(
	dm,
	branching = length(tips) > 1L,
	tips = if (branching) find_tips(dm, root) else root,
	root = random_root(dm)
) {
	check_dpt_suppression(dm)
	
	if (missing(branching) && missing(tips))
		stop('you need to specify at least `branching` or `tips`')
	
	n <- length(dm@phi0)
	stopifnot(is.logical(branching), length(branching) == 1L)
	stopifnot(is.integer(tips), length(tips) %in% c(1L, 3L))
	
	tip_labels <- rep(NA_integer_, n)
	tip_labels[tips] <- tips
	
	if (branching) {
		dpt <- vapply(tips, function(cell) dpt_to_cell(dm, cell), double(n))
		bid <- apply(dpt, 2, order)
		
		# cut it into three branches
		branch <- lapply(seq_len(3), function(b) branchcut(dpt, bid, b))
		
		unassigned <- setdiff(seq_len(3L), Reduce(union, branch, integer()))
		data.frame(
			Branch = organize_branches(branch, unassigned),
			Tips = tip_labels,
			DPT = dpt[, 1],
			DPT = dpt[, 2:3])  # DPT.1 & DPT.2
	} else {
		data.frame(
			Branch = rep(factor('branch 1'), n),
			Tips = tip_labels,
			DPT = dpt_to_cell(dm, tips[[1]]))
	}
}


#' Calculate DPTs
#' 
#' Given a cell index, returns the DPT of this cell to all cells
#' 
#' @param ts    A \code{\link{DiffusionMap}} object
#' @param cell  Index of a cell for which all DPTs should be calculated
#' 
#' @name dpt helpers
#' @export
dpt_to_cell <- function(dm, cell) {
	check_dpt_suppression(dm)
	cell_propagations <- dm@propagations[cell, ]
	
	apply(dm@propagations, 1, function(row)
		sqrt(sum((cell_propagations - row) ^ 2)))
}


check_dpt_suppression <- function(dm) {
	if (is.null(dm@propagations))
		stop('DiffusionMap was created with suppress_dpt = TRUE')
}


# This function organizes the cell arrays branch[[i]] and unassigned (created in dpt)
# to build Branch labels (length <- number of cells) indicating the branch each cell belongs to.
# Cells which are assigned to more than one branch in dpt as well
# as cells which are not assigned to any branch are defined as undeceided (label NA)
#' @importFrom utils combn
organize_branches <- function(branch, unassigned) {
	n <- do.call(max, branch)
	
	intersect_branches <- function(bs) intersect(branch[[bs[[1]]]], branch[[bs[[2]]]])
	branch_intersections <- lapply(combn(3L, 2L, simplify = FALSE), intersect_branches)
	inters <- Reduce(union, branch_intersections, integer())
	
	branch <- lapply(branch, function(b) setdiff(b, inters))
	branch_nums <- seq_along(branch)  # TODO: change
	
	branch_labels <- paste('branch', branch_nums)
	branches_label <- paste(branch_nums, collapse = ',')
	unassigned_label <- paste('unassigned', branches_label)
	uncertain_label <- paste('uncertain', branches_label)
	
	levels <- c(uncertain_label, unassigned_label, branch_labels)
	
	labels <- factor(rep(uncertain_label, n), levels)
	for (b in seq_along(branch)) {
		labels[branch[[b]]] <- paste('branch', branch_nums[[b]])
	}
	labels[unassigned] <- unassigned_label
	labels
}
