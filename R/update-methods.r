#' Update old \link{DiffusionMap}s or \link{Sigmas} to a newer version
#' 
#' @param object  A \link{DiffusionMap} or \link{Sigmas} object created with an older destiny release
#' @param ...     ignored
#' @param verbose tells what is being updated
#' 
#' @return A \link{DiffusionMap} or \link{Sigmas} object that is valid when used with the current destiny release
#' 
#' @aliases updateObject,DiffusionMap-method updateObject,Sigmas-method
#' 
#' @importFrom methods setMethod validObject .hasSlot slot slot<- 
#' @importFrom Matrix Matrix
#' @importFrom BiocGenerics updateObject
#' @name updateObject-methods
#' @export
setMethod('updateObject', 'DiffusionMap', function(object, ..., verbose = FALSE) {
	if (verbose) 
		message("updateObject(object = 'DiffusionMap')")
	
	if (!.hasSlot(object, 'distance'))
		slot(object, 'distance', check = FALSE) <- 'euclidean'
	
	if (!.hasSlot(object, 'transitions'))
		slot(object, 'transitions', check = FALSE) <- Matrix(0, length(object@d), length(object@d), sparse = TRUE)
	
	if (!.hasSlot(object, 'd.norm'))
		slot(object, 'd_norm', check = FALSE) <- rep(NA, length(object@d))
	
	object <- update_slot_names(object, c('data.env', 'd.norm', 'density.norm', 'censor.val', 'censor.range', 'missing.range'))
	
	object@sigmas <- updateObject(object@sigmas)
	
	validObject(object)
	object
})


#' @name updateObject-methods
#' @export
setMethod('updateObject', 'Sigmas', function(object, ..., verbose = FALSE) {
	if (verbose) 
		message("updateObject(object = 'Sigmas')")
	
	object <- update_slot_names(object, c('log.sigmas', 'dim.norms', 'optimal.sigma', 'optimal.idx', 'avrd.norms'))
	
	object
})


update_slot_names <- function(object, old_slots, new_slots = sub('.', '_', old_slots, fixed = TRUE)) {
	if (.hasSlot(object, old_slots[[1]]))
		mapply(function(old_slot, new_slot) {
			slot(object, new_slot, check = FALSE) <- slot(object, old_slot)
			NULL
		}, old_slots, new_slots)
	
	object
}
