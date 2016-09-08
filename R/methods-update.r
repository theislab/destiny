#' @include utils.r
NULL

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
#' @name updateObject-method
#' @export
setMethod('updateObject', 'DiffusionMap', function(object, ..., verbose = FALSE) {
	if (verbose) 
		message("updateObject(object = 'DiffusionMap')")
	
	if (!hasattr(object, 'distance'))
		slot(object, 'distance', check = FALSE) <- 'euclidean'
	
	if (!hasattr(object, 'transitions'))
		slot(object, 'transitions', check = FALSE) <- NULL
	
	if (!hasattr(object, 'd.norm'))  # upgrade only nonexistence, name later
		slot(object, 'd.norm', check = FALSE) <- rep(NA_real_, length(object@d))
	
	if (!hasattr(object, 'n_local'))
		slot(object, 'n_local', check = FALSE) <- 5L
	
	object <- update_slot_names(object, c('data.env', 'd.norm', 'density.norm', 'censor.val', 'censor.range', 'missing.range'))
	
	slot(object, 'sigmas', check = FALSE) <- updateObject(object@sigmas)
	
	validObject(object)
	object
})


#' @name updateObject-method
#' @export
setMethod('updateObject', 'Sigmas', function(object, ..., verbose = FALSE) {
	if (verbose) 
		message("updateObject(object = 'Sigmas')")
	
	object <- update_slot_names(object, c('log.sigmas', 'dim.norms', 'optimal.sigma', 'optimal.idx', 'avrd.norms'))
	
	object
})


update_slot_names <- function(object, old_slots, new_slots = sub('.', '_', old_slots, fixed = TRUE)) {
	atts <- attributes(object)
	update_idx <- old_slots %in% names(atts)
	if (!length(update_idx)) return(object)
	
	slot_idx <- na.omit(match(old_slots[update_idx], names(atts)))
	
	names(atts)[slot_idx] <- new_slots[update_idx]
	attributes(object) <- atts
	
	object
}
