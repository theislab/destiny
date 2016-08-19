#' @include accessor-generics.r
NULL

#' DiffusionMap accession methods
#' 
#' Get and set eigenvalues, eigenvectors, and sigma(s) of a \link{DiffusionMap} object or print information about a DiffusionMap
#' 
#' @param object  A DiffusionMap
#' @param value   Vector of eigenvalues or matrix of eigenvectors to get/set
#' 
#' @return The assigned or retrieved value
#' 
#' @seealso \link{DiffusionMap extraction}, \link{DiffusionMap methods}, \link{coercions} for more methods
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
#' @aliases
#'   eigenvalues,DiffusionMap-method   eigenvectors,DiffusionMap-method   sigmas,DiffusionMap-method
#' eigenvalues<-,DiffusionMap-method eigenvectors<-,DiffusionMap-method sigmas<-,DiffusionMap-method
#'   dataset,DiffusionMap-method   distance,DiffusionMap-method   optimal_sigma,DiffusionMap-method
#' dataset<-,DiffusionMap-method distance<-,DiffusionMap-method
#' 
#' @importFrom methods is setGeneric
#' @name DiffusionMap accessors
#' @include sigmas.r
#' @include diffusionmap.r
NULL


#' @name DiffusionMap accessors
#' @export
setMethod('eigenvalues', 'DiffusionMap', function(object) object@eigenvalues)

#' @name DiffusionMap accessors
#' @export
setMethod('eigenvalues<-', 'DiffusionMap', function(object, value) {
	object@eigenvalues <- value
	object
})


#' @name DiffusionMap accessors
#' @export
setMethod('eigenvectors', 'DiffusionMap', function(object) object@eigenvectors)

#' @name DiffusionMap accessors
#' @export
setMethod('eigenvectors<-', 'DiffusionMap', function(object, value) {
	object@eigenvectors <- value
	validObject(object)
	object
})


#' @name DiffusionMap accessors
#' @export
setMethod('sigmas', 'DiffusionMap', function(object) object@sigmas)

#' @name DiffusionMap accessors
#' @export
setMethod('sigmas<-', 'DiffusionMap', function(object, value) {
	object@sigmas <- value
	validObject(object)
	object
})


#' @name DiffusionMap accessors
#' @export
setMethod('dataset', 'DiffusionMap', function(object) object@data_env$data)

#' @name DiffusionMap accessors
#' @export
setMethod('dataset<-', 'DiffusionMap', function(object, value) {
	object@data_env$data <- value
	validObject(object)
	object
})


#' @name DiffusionMap accessors
#' @export
setMethod('distance', 'DiffusionMap', function(object) object@distance)

#' @name DiffusionMap accessors
#' @export
setMethod('distance<-', 'DiffusionMap', function(object, value) {
	object@distance <- value
	validObject(object)
	object
})


#' @name DiffusionMap accessors
#' @export
setMethod('optimal_sigma', 'DiffusionMap', function(object) optimal_sigma(sigmas(object)))
