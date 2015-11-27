#' @importFrom BiocGenerics updateObject
NULL

#' Update old DiffusionMaps to a newer version
#' 
#' @param object  A DiffusionMap/Sigmas object
#' @param ...     ignored
#' @param verbose tells what is being updated
#' 
#' @aliases updateObject,Diffusionmap-method updateObject,sigmas-method
#' @name destiny updateObject methods
NULL

#' @name destiny updateObject methods
#' @export
setMethod('updateObject', 'DiffusionMap', function(object, ..., verbose = FALSE) {
	if (verbose) 
		message("updateObject(object = 'DiffusionMap')")
	
	if (!.hasSlot(object, 'distance'))
		slot(object, 'distance', check = FALSE) <- 'euclidean'
	
	validObject(object)
	object
})
