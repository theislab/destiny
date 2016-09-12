#' @include dpt.r
NULL

#' DPT methods
#' 
#' Methods for the \link{DPT} class. \code{branch_divide} subdivides branches for plotting (see the examples).
#' 
#' @param dpt,object  DPT object
#' @param divide      Vector of branch numbers to use for division
#' @param value       Value of slot to set
#' 
#' @return \code{branch_divide} and \code{dataset<-} return the changed object, \code{dataset} the extracted data.
#' 
#' @examples
#' data(guo_norm)
#' dpt <- DPT(DiffusionMap(guo_norm))
#' dpt_9_branches <- branch_divide(dpt, 1:3)
#' plot(dpt_9_branches, col_by = 'branch')
#' 
#' @seealso \link{plot.DPT} uses \code{branch_divide} for its \code{divide} argument.
#' 
#' @aliases dataset,DPT-method dataset<-,DPT-method
#' @importFrom stats na.omit
#' @name DPT methods
#' @export
branch_divide <- function(dpt, divide = integer(0L)) {
	if (!is(dpt, 'DPT')) stop('branch_divide needs to be called on a DPT object, not a ', class(dpt))
	if (length(divide) == 0L) return(dpt)
	
	for (b in divide) {
		super_rows <- dpt@branch[, 1] == b & !is.na(dpt@branch[, 1])
		if (!any(super_rows)) {
			available <- na.omit(unique(dpt@branch[, 1]))
			stop('invalid branch to divide ', b, ' not in ', available)
		}
		
		# shift sub branches/tips to the left
		dpt@branch[super_rows, ] <- cbind(dpt@branch[super_rows, -1], NA)
		dpt@tips  [super_rows, ] <- cbind(dpt@tips  [super_rows, -1], NA)
		
		# TODO: maybe also modify DPT?
	}
	
	vacant_levels <- apply(dpt@branch, 2L, function(col) all(is.na(col)))
	dpt@branch <- dpt@branch[, !vacant_levels]
	dpt@tips   <- dpt@tips  [, !vacant_levels]
	
	dpt
}

#' @name DPT methods
#' @export
setMethod('dataset', 'DPT', function(object) dataset(object@dm))

#' @name DPT methods
#' @export
setMethod('dataset<-', 'DPT', function(object, value) {
	dataset(object@dm) <- value
	validObject(object)
	object
})