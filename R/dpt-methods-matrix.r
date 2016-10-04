#' @include dpt.r s4-unions.r
NULL

#' DPT Matrix methods
#' 
#' Treat DPT object as a matrix of cell-by-cell DPT distances.
#' 
#' @aliases [.DPT nrow.DPT ncol.DPT dim.DPT
#' @name DPT matrix methods
#' @export
setMethod('[', c('DPT', 'index', 'index', 'logicalOrMissing'), function(x, i, j, ..., drop = TRUE) {
	evas <- eigenvalues(x@dm)
	eves <- eigenvectors(x@dm)
	
	# get numeric from negative or logical indices
	i <- seq_len(nrow(x))[i]
	j <- seq_len(nrow(x))[j]
	
	norm <- array(NA, c(length(i), length(j), length(evas)))
	ev_dist <- norm
	
	for (e in seq_along(evas)) norm[, , e] <- evas[[e]]
	norm <- norm / (1-norm)
	
	do.call(mapply, c(list(function(ii, jj) {
		ev_dist[ii, jj, ] <<- eves[i[[ii]], ] - eves[j[[jj]], ]
	}), expand.grid(ii = seq_along(i), jj = seq_along(j))))
	
	r <- sqrt(apply(norm^2 * ev_dist^2, 1:2, sum))
	if (drop && 1L %in% dim(r)) dim(r) <- NULL
	r
})

#' @name DPT matrix methods
#' @export
setMethod('[', c('DPT', 'index', 'missing', 'logicalOrMissing'), function(x, i, j, ..., drop = TRUE) {
	x[i, seq_len(nrow(x)), ..., drop = drop]
})

#' @name DPT matrix methods
#' @export
setMethod('[', c('DPT', 'missing', 'index', 'logicalOrMissing'), function(x, i, j, ..., drop = TRUE) {
	x[seq_len(nrow(x)), j, ..., drop = drop]
})

#' @name DPT matrix methods
#' @export
setMethod('[', c('DPT', 'missing', 'missing', 'logicalOrMissing'), function(x, i, j, ..., drop = TRUE) {
	x[seq_len(nrow(x)), seq_len(nrow(x)), ..., drop = drop]
})

#' @name DPT matrix methods
#' @export
setMethod('nrow', 'DPT', function(x) length(x@dm@d))
#' @name DPT matrix methods
#' @export
setMethod('ncol', 'DPT', function(x) length(x@dm@d))
#' @name DPT matrix methods
#' @export
setMethod('dim', 'DPT', function(x) c(nrow(x), ncol(x)))
