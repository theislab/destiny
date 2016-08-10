#' DiffusionMap methods
#' 
#' Get and set eigenvalues, eigenvectors, and sigma(s) of a \link{DiffusionMap} object or print information about a DiffusionMap
#' 
#' @param dm,x,object  A DiffusionMap
#' @param value  Vector of eigenvalues or matrix of eigenvectors to get/set
#' 
#' @return The assigned or retrieved value
#' 
#' @seealso \code{\link{DiffusionMap}}, \code{\link{DiffusionMap extraction}}
#' 
#' @examples
#' data(guo)
#' dm <- DiffusionMap(guo)
#' eigenvalues(dm)
#' eigenvectors(dm)
#' sigmas(dm)
#' dataset(dm)
#' optimal_sigma(dm)
#' 
#' @aliases DiffusionMap-methods values values<- vectors vectors<- sigmas sigmas<- optimal_sigma,DiffusionMap-method print,DiffusionMap-method show,DiffusionMap-method
#' 
#' @importFrom methods is
#' @name DiffusionMap methods
#' @include sigmas.r
#' @include diffusionmap.r
NULL


#' @name DiffusionMap methods
#' @export
eigenvalues <- function(dm) {
	stopifnot(is(dm, 'DiffusionMap'))
	dm@eigenvalues
}

#' @name DiffusionMap methods
#' @export
`eigenvalues<-` <- function(dm, value) {
	stopifnot(is(dm, 'DiffusionMap'))
	dm@eigenvalues <- value
	validObject(dm)
	dm
}


#' @name DiffusionMap methods
#' @export
eigenvectors <- function(dm) {
	stopifnot(is(dm, 'DiffusionMap'))
	dm@eigenvectors
}

#' @name DiffusionMap methods
#' @export
`eigenvectors<-` <- function(dm, value) {
	stopifnot(is(dm, 'DiffusionMap'))
	dm@eigenvectors <- value
	validObject(dm)
	dm
}


#' @name DiffusionMap methods
#' @export
sigmas <- function(dm) {
	stopifnot(is(dm, 'DiffusionMap'))
	dm@sigmas
}

#' @name DiffusionMap methods
#' @export
`sigmas<-` <- function(dm, value) {
	stopifnot(is(dm, 'DiffusionMap'))
	dm@sigmas <- value
	validObject(dm)
	dm
}


#' @name DiffusionMap methods
#' @export
dataset <- function(dm) {
	stopifnot(is(dm, 'DiffusionMap'))
	dm@data_env$data
}

#' @name DiffusionMap methods
#' @export
`dataset<-` <- function(dm, value) {
	stopifnot(is(dm, 'DiffusionMap'))
	dm@data_env$data <- value
	validObject(dm)
	dm
}


#' @name DiffusionMap methods
#' @export
setMethod('optimal_sigma', 'DiffusionMap', function(object) sigmas(object)@optimal_sigma)


#' @importFrom utils str
#' 
#' @name DiffusionMap methods
#' @export
setMethod('print', 'DiffusionMap', function(x) {
	cat(sprintf('DiffusionMap (%s Diffusion components and %s samples)\n', length(eigenvalues(x)), nrow(eigenvectors(x))))
	cat('eigenvalues:  '); str(eigenvalues(x))
	cat('eigenvectors: '); str(structure(eigenvectors(x), dimnames = NULL))
	cat('  ..colnames: '); str(colnames(eigenvectors(x)), vec.len = 4)
	cat(sprintf('optimal_sigma: %s\n', optimal_sigma(x)))
	cat(sprintf('@distance:     %s', x@distance));
	invisible(x)
})

#' @name DiffusionMap methods
#' @export
setMethod('show', 'DiffusionMap', function(object) {
	print(object)
	invisible()
})
