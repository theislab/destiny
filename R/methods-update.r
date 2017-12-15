#' @include utils.r
NULL

#' Update old destiny objects to a newer version.
#' 
#' Handles \link{DiffusionMap}, \link{Sigmas}, and \link{GeneRelevance}.
#' 
#' @param object  An object created with an older destiny release
#' @param ...     ignored
#' @param verbose tells what is being updated
#' 
#' @return A \link{DiffusionMap} or \link{Sigmas} object that is valid when used with the current destiny release
#' 
#' @aliases updateObject,DiffusionMap-method updateObject,Sigmas-method updateObject,GeneRelevance-method
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
	
	if (!hasattr(object, 'rotate'))
		slot(object, 'rotate', check = FALSE) <- TRUE  # old ones were rotated by default
	
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


#' @name updateObject-method
#' @export
setMethod('updateObject', 'GeneRelevance', function(object, ..., verbose = FALSE) {
	if (verbose) 
		message("updateObject(object = 'GeneRelevance')")
	
	if (!hasattr(object, 'distance'))
		slot(object, 'distance', check = FALSE) <- 'euclidean'
	
	# the dimensions were switched to fit the convention elsewhere in the package.
	if (!hasattr(object, 'smooth_window')) {
		object@partials <- aperm(object@partials, c(2, 1, 3))
		object@partials_norm <- t(object@partials_norm)
		slot(object, 'smooth_window', check = FALSE) <- NA_real_
		slot(object, 'smooth_alpha', check = FALSE) <- NA_real_
	}
	
	object
})


#' @importFrom stats na.omit
update_slot_names <- function(object, old_slots, new_slots = sub('.', '_', old_slots, fixed = TRUE)) {
	atts <- attributes(object)
	update_idx <- old_slots %in% names(atts)
	if (!length(update_idx)) return(object)
	
	slot_idx <- na.omit(match(old_slots[update_idx], names(atts)))
	
	names(atts)[slot_idx] <- new_slots[update_idx]
	attributes(object) <- atts
	
	object
}
