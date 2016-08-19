#' DiffusionMap extraction methods
#' 
#' @param x       \code{\link{DiffusionMap}} object
#' @param i,name  Name of a diffusion component (DC1,...), or column from the data
#' @param j       N/A
#' @param ...     ignored
#' 
#' @seealso \link[base]{Extract}, \code{\link[base]{names}} for the generics. \link{DiffusionMap accessors}, \link{DiffusionMap methods}, \link{coercions} for more methods
#' 
#' @name DiffusionMap extraction
#' @aliases names.DiffusionMap names,DiffusionMap-method $,DiffusionMap-method [[,DiffusionMap,character,missing-method
NULL


#' @importFrom methods is
#' @importFrom Biobase featureNames varLabels
#' @name DiffusionMap extraction
#' @export
setMethod('names', 'DiffusionMap', function(x) {
	dta <- dataset(x)
	data_names <- if (is(dta, 'ExpressionSet')) {
		c(featureNames(dta), varLabels(dta))
	} else {
		colnames(dta)
	}
	c(colnames(eigenvectors(x)), data_names)
})


#' @importFrom methods is
#' @importFrom Biobase exprs featureNames
#' @name DiffusionMap extraction
#' @export
setMethod('[[', c('DiffusionMap', 'character', 'missing'), function(x, i, j, ...) {
	dta <- if (grepl('^DC\\d+$', i)) eigenvectors(x) else dataset(x)
	if (is.matrix(dta)) {
		dta[, i]
	} else if (is(dta, 'ExpressionSet') && i %in% featureNames(dta)) {
		exprs(dta)[i, ]
	} else {  # data.frame or phenoData
		dta[[i]]
	}
})


#' @name DiffusionMap extraction
#' @export
setMethod('$', 'DiffusionMap', function(x, name) x[[name]])
