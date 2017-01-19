#' Projection distance
#' 
#' @param dists_or_data  If new_data is not preset, this needs to be a \code{\link[=proxy::dist]{crossdist}} object
#' @param new_data       A \code{\link[base]{matrix}} \code{\link[Matrix]{Matrix}}, or \code{\link[base]{data.frame}}.
#'                       If given, \code{dists_or_data e} needs to be of one of the same classes.
#' @param ...            Passed to \code{\link{proxa::dist}}
#' 
#' @export
setGeneric(
	'projection_dist',
	function(dists_or_data, new_data, ...) standardGeneric('projection_dist'),
	signature = c('dists_or_data', 'new_data'))

.projection_dist_impl <- function(dists) {
	
}

#' @export
setMethod('projection_dist', c('anyMatrix', 'anyMatrix'), function(dists_or_data, new_data, ...)
	.projection_dist_impl(proxy::dist(dists_or_data, new_data, ...)))

#' @export
setMethod('projection_dist', c('crossdist', 'missing'), function(dists_or_data, new_data, ...)
	.projection_dist_impl(dists_or_data))
