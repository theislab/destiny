#' @include diffusionmap.r
NULL

#' @importFrom Hmisc rcorr
rankcor.dist <- function(x, y) {
	if (length(x) != length(y))
		stop('Error: The vectors have different length!')
	
	1 - rcorr(x, y, type = 'spearman')$r[1, 2]
}

centered.cosine.dist <- function(x, y) 1 - cor(x, y)

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
	#censor imples euclidean distance
	
	sigma <- optimal.sigma(dm)
	if (!censor) {
		measure <- switch(dm@distance, euclidean = 'Euclidean', cosine = centered.cosine.dist, rankcor = rankcor.dist, stop('Unknown distance measure'))
		d2 <- unclass(proxy::dist(new.data, data, measure) ^ 2) # matrix (dense)
		
		#TODO: zeros not on diag
		trans.p <- exp(-d2 / (2 * sigma ^ 2))
		trans.p[d2 == 0] <- 0
	} else {
		trans.p <- predict.censoring(data, new.data, dm@censor.val, dm@censor.range, dm@missing.range, sigma)
	}
	
	#trans.p:  columns: old data, rows: new data
	d_new <- rowSums(trans.p, na.rm = TRUE)
	norm_p <- get_norm_p(trans.p, dm@d, d_new, dm@density.norm)
	rm(trans.p)  # free memory
	
	d_norm_new <- rowSums(norm_p)
	
	#
	Hp_new <- rotate_norm_p(norm_p, d_norm_new)
	rm(norm_p)  # free memory
	
	phi <- cbind(dm@eigenvec0, eigenvectors(dm))
	
	eig.vec.norm <- Hp_new %*% phi %*% Diagonal(x = c(1, 1 / eigenvalues(dm)))
	eig.vec.norm[, -1]
	
	#new version doesn’t quite work yet
	# d_rot     <- Diagonal(x = dm@d_norm  ^ -.5)
	# d_rot_new <- Diagonal(x = d_norm_new ^ -.5)
	# M_new <- d_rot_new %*% norm_p %*% d_rot
	# 
	# t(t(M_new %*% eigenvectors(dm)) / eigenvalues(dm))
}
