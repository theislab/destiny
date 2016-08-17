#' Compatibility versions of functions for destiny 1.x scripts
#' 
#' This is most useful if you have old scripts that use the old argument order and the default being global sigma.
#' Simply put \code{DiffusionMap <- DiffusionMapCompat} on top of your script to make it work with destiny 2.x.
#' 
#' @param data,n,n.eigs,start,censor.val,censor.range,missing.range,vars,verbose,... Corresponding parameters in multiple functions
#' @param sigma,k,density.norm,distance                                              Corresponding \code{\link{DiffusionMap}} parameters
#' @param r,hue,gamma,light,dark,reverse                                             Corresponding \code{\link{cube_helix}} parameters
#' @param dm,new.data                                                                Corresponding \code{\link{dm_predict}} parameters
#' @param M,sym                                                                      Corresponding \code{\link{eig_decomp}} parameters
#' @param min.k,small,big                                                            Corresponding \code{\link{find_dm_k}} parameters
#' @param step.size,steps,sample.rows,early.exit                                     Corresponding \code{\link{find_sigmas}} parameters
#' @param object                                                                     Corresponding \code{\link{optimal_sigma}} parameters
#' 
#' @name destiny 1 compatibility
#' @importFrom BiocGenerics updateObject
#' @keywords internal
#' @export
DiffusionMapCompat <- function(
	data,
	sigma = 'global',
	k = find_dm_k(nrow(data) - 1L),
	n.eigs = min(20L, nrow(data) - 2L),
	density.norm = TRUE,
	...,
	distance = c('euclidean', 'cosine', 'rankcor'),
	censor.val = NULL, censor.range = NULL,
	missing.range = NULL,
	vars = NULL,
	verbose = !is.null(censor.range)
) DiffusionMap(
	data = data,
	sigma = sigma,
	k = k,
	n_eigs = n.eigs,
	density_norm = density.norm,
	...,
	distance = distance,
	censor_val = censor.val, censor_range = censor.range,
	missing_range = missing.range,
	vars = vars,
	verbose = verbose
)

#' @name destiny 1 compatibility
#' @export
cube.helix <- function(n = 6, start = 0, r = .4, hue = .8, gamma = 1, light = .85, dark = .15, reverse = FALSE) {
	.Deprecated('cube_helix')
	cube_helix(n, start, r, hue, gamma, light, dark, reverse)
}

#' @name destiny 1 compatibility
#' @export
dm.predict <- function(dm, new.data, verbose = FALSE) {
	.Deprecated('dm_predict')
	dm_predict(updateObject(dm), new.data, verbose)
}

#' @name destiny 1 compatibility
#' @export
eig.decomp <- function(M, n.eigs, sym = isSymmetric(M)) {
	.Deprecated('eig_decomp')
	eig_decomp(M, n.eigs, sym)
}

#' @name destiny 1 compatibility
#' @export
find.dm.k <- function(n, min.k = 100L, small = 1000L, big = 10000L) {
	.Deprecated('find_dm_k')
	find_dm_k(n, min.k, small, big)
}

#' @name destiny 1 compatibility
#' @export
find.sigmas <- function(
	data,
	step.size = 0.1,
	steps = 10L,
	start = NULL,
	sample.rows = 500L,
	early.exit = FALSE,
	...,
	censor.val = NULL, censor.range = NULL,
	missing.range = NULL,
	vars = NULL,
	verbose = TRUE
) {
	.Deprecated('find_sigmas')
	find_sigmas(
		data = data,
		step_size = step.size, steps = steps,
		start = start,
		sample_rows = sample.rows,
		early_exit = early.exit,
		...,
		censor_val = censor.val, censor_range = censor.range,
		missing_range = missing.range,
		vars = vars,
		verbose = verbose
	)
}

#' @name destiny 1 compatibility
#' @export
optimal.sigma <- function(object) {
	.Deprecated('optimal_sigma')
	optimal_sigma(updateObject(object))
}
