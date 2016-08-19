#' Diffusion Pseudo Time
#'
#' Create pseudotime ordering and assigns cell to one of three branches
#' 
#' @param dm           A \code{\link{DiffusionMap}} object. Its transition probabilities will be used to calculate the DPT
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
#' @slot branch  Branch labels for each cell; \code{1:3} or \code{NA} for undeceided
#' @slot parent  \code{\link{matrix}} of parent branches (may be of \code{ncol(...) == 0})
#' @slot tips    Indices of tips
#' @slot dpt     Diffusion pseudotime in respect to the root cell (and other tips if \code{branching == TRUE})
#' @slot dm      \code{\link{DiffusionMap}} used to create this DPT object
#' 
#' @aliases DPT-class
#' 
#' @importFrom methods setClass
#' @name DPT
#' @export
setClass(
	'DPT',
	slots = c(
		branch = 'integer',
		parent = 'matrix', # 'integer'
		tips   = 'integer',
		dpt    = 'matrix', # 'double'
		dm     = 'DiffusionMap'),
	validity = function(object) {
		TRUE
	})

#' @name DPT
#' @export
DPT <- function(
	dm,
	branching = length(tips) > 1L,
	tips = if (branching) find_tips(dm, root) else root,
	root = random_root(dm)
) {
	if (!is(dm, 'DiffusionMap')) stop('dm needs to be of class DiffusionMap, not ', class(dm))
	if (missing(branching) && missing(tips))
		stop('you need to specify at least `branching` or `tips`')
	
	n <- length(dm@d_norm)
	stopifnot(is.logical(branching), length(branching) == 1L)
	stopifnot(is.integer(tips), length(tips) %in% c(1L, 3L))
	
	if (branching) {
		dpt <- vapply(tips, function(cell) dpt_to_cell(dm, cell), double(n))
		bid <- apply(dpt, 2, order)
		
		# cut it into three branches
		branch <- lapply(seq_len(3), function(b) branchcut(dpt, bid, b))
		
		new(
			'DPT',
			branch = organize_branches(branch),
			parent = matrix(NA_integer_, n, 0L),
			tips = tips,
			dpt = dpt,
			dm = dm)
	} else {
		new(
			'DPT',
			branch = rep(1L, n),
			parent = matrix(NA_integer_, n, 0L),
			tips = tips,
			dpt = dpt_to_cell(dm, tips[[1]]),
			dm = dm)
	}
}


#' Calculate DPTs
#' 
#' Given a cell index, returns the DPT of this cell to all cells
#' 
#' @param dm    A \code{\link{DiffusionMap}} object
#' @param cell  Index of a cell for which all DPTs should be calculated
#' 
#' @name dpt helpers
#' @export
dpt_to_cell <- function(dm, cell) {
	propagations <- propagation_matrix(dm)
	cell_propagations <- propagations[cell, ]
	
	apply(propagations, 1, function(row)
		sqrt(sum((cell_propagations - row) ^ 2)))
}
