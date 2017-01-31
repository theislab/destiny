#' Projection distance
#' 
#' @param dm        A \code{\link{DiffusionMap}} object.
#' @param new_data  New data points to project into the diffusion map.
#'                  Can be a \link[base]{matrix}, \link[base]{data.frame}, or an \link[Biobase]{ExpressionSet}.
#' @param ...       Passed to \code{\link{proxy::dist}}
#' 
#' @importFrom FNN get.knnx
#' @export
projection_dist <- function(dm, new_data, ..., verbose = FALSE) {
	predicted <- dm_predict(dm, new_data, ..., verbose = verbose)
	
	evs <- eigenvectors(dm)
	
	nns <- get.knnx(evs, predicted, 1)
	
	#nn_idx <- nns$nn.index[, 1]
	nn_dist <- nns$nn.dist[, 1]
	
	nn_dist
}
