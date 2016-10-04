#' Diffusion Pseudo Time
#'
#' Create pseudotime ordering and assigns cell to one of three branches
#' 
#' @param dm       A \code{\link{DiffusionMap}} object. Its transition probabilities will be used to calculate the DPT
#' @param tips     The cell index/indices from which to calculate the DPT(s) (integer of length 1-3)
#' @param ...      All parameters after this have to be specified by name
#' @param w_width  Window width to use for deciding the branch cutoff
#' 
#' @return A \code{DPT} object:
#' 
#' @slot branch  Branch labels for each cell; \code{1:3} or \code{NA} for undeceided
#' @slot tips    \code{\link[base]{matrix}} indicating if a cell is a tip of the corresponding banch level
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
		dm     = 'DiffusionMap'),
	validity = function(object) {
		TRUE  # TODO
	})

#' @name DPT
#' @export
DPT <- function(dm, tips = random_root(dm), ..., w_width = .1) {
	if (!is(dm, 'DiffusionMap')) stop('dm needs to be of class DiffusionMap, not ', class(dm))
	if (!length(tips) %in% 1:3) stop('you need to specify 1-3 tips, got ', length(tips))
	
	dpt <- dummy_dpt(dm)
	
	stats <- tipstats(dpt, tips)
	branches <- auto_branch(dpt, seq_len(nrow(dpt)), stats, w_width)
	
	colnames(branches$branch) <- paste0('Branch', seq_len(ncol(branch)))
	colnames(branches$tips)   <- paste0('Tips',   seq_len(ncol(tip_mat)))
	
	dpt@branch <- branches$branch
	dpt@tips <- branches$tips
	dpt
}

dummy_dpt <- function(dm) new('DPT', branch = matrix(), tips = matrix(), dm = dm)
