#' @include accessor-generics.r
NULL

#' DiffusionMap accession methods
#' 
#' Get and set eigenvalues, eigenvectors, and sigma(s) of a \link{DiffusionMap} object.
#' 
#' @param object  A DiffusionMap
#' @param value   Vector of eigenvalues or matrix of eigenvectors to get/set
#' 
#' @return The assigned or retrieved value
#' 
#' @seealso \link{extractions}, \link{DiffusionMap methods}, \link{coercions} for more methods
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
#' @importFrom methods is setGeneric
#' @name DiffusionMap accession methods
#' @rdname DiffusionMap-accessors
#' @include sigmas.r
#' @include diffusionmap.r
NULL


#' @rdname DiffusionMap-accessors
#' @export
setMethod('eigenvalues', 'DiffusionMap', function(object) object@eigenvalues)

#' @rdname DiffusionMap-accessors
#' @export
setMethod('eigenvalues<-', 'DiffusionMap', function(object, value) {
	object@eigenvalues <- value
	object
})


#' @rdname DiffusionMap-accessors
#' @export
setMethod('eigenvectors', 'DiffusionMap', function(object) object@eigenvectors)

#' @rdname DiffusionMap-accessors
#' @export
setMethod('eigenvectors<-', 'DiffusionMap', function(object, value) {
	object@eigenvectors <- value
	validObject(object)
	object
})


#' @rdname DiffusionMap-accessors
#' @export
setMethod('sigmas', 'DiffusionMap', function(object) object@sigmas)

#' @rdname DiffusionMap-accessors
#' @export
setMethod('sigmas<-', 'DiffusionMap', function(object, value) {
	object@sigmas <- value
	validObject(object)
	object
})


#' @rdname DiffusionMap-accessors
#' @export
setMethod('dataset', 'DiffusionMap', function(object) object@data_env$data)

#' @rdname DiffusionMap-accessors
#' @export
setMethod('dataset<-', 'DiffusionMap', function(object, value) {
	object@data_env$data <- value
	validObject(object)
	object
})


#' @rdname DiffusionMap-accessors
#' @export
setMethod('distance', 'DiffusionMap', function(object) object@distance)

#' @rdname DiffusionMap-accessors
#' @export
setMethod('distance<-', 'DiffusionMap', function(object, value) {
	object@distance <- value
	validObject(object)
	object
})


#' @rdname DiffusionMap-accessors
#' @export
setMethod('optimal_sigma', 'DiffusionMap', function(object) optimal_sigma(sigmas(object)))
