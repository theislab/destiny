#' Projection distance
#' 
#' @param dm        A \code{\link{DiffusionMap}} object.
#' @param new_dcs   Diffusion component matrix of which to calculate the distance to the data.
#' @param ...       Passed to \code{\link{proxy::dist}} if \code{new_data} was passed.
#' @param new_data  New data points to project into the diffusion map.
#'                  Can be a \link[base]{matrix}, \link[base]{data.frame}, or an \link[Biobase]{ExpressionSet}.
#' 
#' @importFrom FNN get.knnx
#' @export
projection_dist <- function(dm, new_dcs = NULL, ..., new_data, verbose = FALSE) {
	if (is.null(new_dcs))
		new_dcs <- dm_predict(dm, new_data, ..., verbose = verbose)
	else if (!missing(new_data))
		stop('only pass one of new_dcs and new_data')
	
	evs <- eigenvectors(dm)
	
	nns <- get.knnx(evs, new_dcs, 1)
	
	#nn_idx <- nns$nn.index[, 1]
	nn_dist <- nns$nn.dist[, 1]
	
	nn_dist
}
