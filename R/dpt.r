#' Diffusion Pseudo Time
#'
#' Create pseudotime ordering and assigns cell to one of three branches
#' 
#' @param dm           A \code{\link{DiffusionMap}} object. Its transition probabilities will be used to calculate the DPT
#' @param branching    Detect a branching? (\code{TRUE} or \code{FALSE})
#' @param root         The root index from which to calculate the DPTs (integer of length 1)
#' 
#' @slot branch  Branch labels for each cell; \code{1:3} or \code{NA} for undeceided
#' @slot parent  \code{\link{matrix}} of parent branches (may be of \code{ncol(...) == 0})
#' @slot dpt     Diffusion pseudotime in respect to the root cell (and other tips if \code{branching == TRUE})
#' @slot dm      \code{\link{DiffusionMap}} used to create this DPT object
#' 
#' @aliases DPT-class
#' 
#' @examples
#' data(guo_norm)
#' dm <- DiffusionMap(guo_norm)
#' dpt <- DPT(dm, branching = TRUE)
#' str(dpt)
#' 
#' @importFrom methods setClass
#' @name DPT
#' @export
setClass(
	'DPT',
	slots = c(
		branch = 'matrix', # 'integer'
		tips   = 'matrix', # 'logical'
		dpt    = 'matrix', # 'double'
		dm     = 'DiffusionMap'),
	validity = function(object) {
		TRUE  # TODO
	})

#' @name DPT
#' @export
DPT <- function(dm, branching = TRUE, root = random_root(dm)) {
	if (!is(dm, 'DiffusionMap')) stop('dm needs to be of class DiffusionMap, not ', class(dm))
	
	n <- length(dm@d_norm)
	stopifnot(is.logical(branching), length(branching) == 1L)
	
	if (branching) {
		propagations <- propagation_matrix(dm)
		stats <- tipstats(propagations, root)
		branches <- auto_branch(propagations, stats)
		
		branch <- branches$branch
		tips <- branches$tips
		dpt <- branches$dpt
	} else {
		propagations <- propagation_matrix(dm)
		dpt_to_root <- dpt_to_cell(propagations, root)
		
		branch <- matrix(rep(1L, n), 1L, n)
		tips <- root
		dpt <- matrix(dpt_to_root, n, 1L)
	}
	
	colnames(branch) <- paste0('Branch', seq_len(ncol(branch)))
	colnames(tips)   <- paste0('Tips',   seq_len(ncol(tips)))
	colnames(dpt)    <- paste0('DPT',    seq_len(ncol(dpt)))
	
	new('DPT', branch = branch, tips = tips, dpt = dpt, dm = dm)
}


dpt_to_cell <- function(propagations, cell) {
	cell_propagations <- propagations[cell, ]
	
	apply(propagations, 1, function(row)
		sqrt(sum((cell_propagations - row) ^ 2)))
}
