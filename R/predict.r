#' @include diffusionmap.r
NULL

#' Predict new data points using an existing DiffusionMap. The resulting matrix can be used in \link[=plot.DiffusionMap]{the plot method for the DiffusionMap}
#' 
#' @param dm        A \link{DiffusionMap} object
#' @param new.data  New data points to project into the diffusion map. Can be a \link[base]{matrix}, \link[base]{data.frame}, or an \link[Biobase]{ExpressionSet}.
#' 
#' @return A \eqn{nrow(new.data) \times ncol(eigenvectors(dif))} matrix of projected diffusion components for the new data.
#' 
#' @examples
#' data(guo)
#' g1 <- guo[, guo$num.cells != 32L]
#' g2 <- guo[, guo$num.cells == 32L]
#' dm <- DiffusionMap(g1)
#' dc2 <- dm.predict(dm, g2)
#' plot(dm, new.dcs = dc2)
#' 
#' @importFrom proxy dist
#' @export
dm.predict <- function(dm, new.data) {
	if (!is(dm, 'DiffusionMap')) stop('dm needs to be a DiffusionMap')
	
	data <- extract.doublematrix(dataset(dm), dm@vars)
	new.data <- extract.doublematrix(new.data, dm@vars)
	if (!ncol(data) == ncol(new.data)) stop('new data needs to have the same features as the one used to create the diffusion map')
	
	censor <- !is.null(dm@censor.val) || any(is.na(data)) || any(is.na(new.data))
	sigma <- optimal.sigma(dm)
	if (!censor) {
		d2 <- unclass(proxy::dist(new.data, data) ^ 2) # matrix (dense)
		
		#TODO: zeros not on diag
		trans.p <- exp(-d2 / (2 * sigma ^ 2))
		trans.p[d2 == 0] <- 0
	} else {
		trans.p <- predict.censoring(data, new.data, dm@censor.val, dm@censor.range, dm@missing.range, sigma)
	}
	
	#columns: old data, rows: new data
	
	trans.p <- as(trans.p, 'dgTMatrix')
	
	max.dist <- max(trans.p@x, na.rm = TRUE)
	stopifsmall(max.dist)
	
	d.new <- rowSums(trans.p, na.rm = TRUE)
	
	if (dm@density.norm) {
		# faster version of H <- trans.p / outer(d.new, dm@d)
		H <- sparseMatrix(trans.p@i, trans.p@j, x = trans.p@x / (d.new[trans.p@i + 1] * dm@d[trans.p@j + 1]), dims = dim(trans.p), index1 = FALSE)
		#creates a dgCMatrix
	} else {
		H <- trans.p
	}
	rm(trans.p)  # free memory
	
	# calculate the inverse of a diagonal matrix by inverting the diagonal
	D <- Diagonal(x = rowSums(H) ^ -1)
	
	Hp <- D %*% H
	rm(H)  # free memory
	
	phi <- cbind(dm@eigenvec0, eigenvectors(dm))
	
	eig.vec.norm <- Hp %*% phi %*% Diagonal(x = c(1, 1 / eigenvalues(dm)))
	eig.vec.norm[, -1]
}
