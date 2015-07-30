#' Logical which
#' 
#' Inverse of \link[base]{which}. Converts an array of numeric or character indices to a logical index array.
#' This function is useful if you need to perform logical operation on an index array but are only given numeric indices.
#' 
#' Either \code{nms} or \code{len} has to be specified.
#' 
#' @param idx       Numeric or character indices.
#' @param nms       Array of names or a sequence. Required if \code{idx} is a character array
#' @param len       Length of output array. Alternative to \code{nms} if \code{idx} is numeric
#' @param useNames  Use the names of nms or idx
#' 
#' @return Logical vector of length \code{len} or the same length as \code{nms}
#' 
#' @examples
#' all(lWhich(2, len = 3) == c(FALSE, TRUE, FALSE))
#' all(lWhich(c('a', 'c'), letters[1:3]) == c(TRUE, FALSE, TRUE))
#' 
#' @export
lWhich <- function(idx, nms = seq_len(len), len = length(nms), useNames = TRUE) {
	rv <- logical(len)
	if (is.character(nms)) # we need names here so that rv[idx] works
		names(rv) <- nms
	
	if (useNames && !is.null(names(idx)))
		names(rv)[idx] <- names(idx)
	
	rv[idx] <- TRUE
	
	if (!useNames) # if we don't want names, we'll remove them if we added them before
		names(rv) <- NULL
	rv
}
