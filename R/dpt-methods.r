#' DPT methods
#' 
#' Methods for the \link{DPT} class
#' 
#' @name DPT methods
NULL

#' @name DPT methods
#' @export
setMethod('dataset', 'DPT', function(object) dataset(object@dm))

#' @name DPT methods
#' @export
setMethod('dataset<-', 'DPT', function(object, value) {
	dataset(object@dm) <- value
	validObject(object)
	object
})