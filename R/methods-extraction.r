#' Extraction methods
#' 
#' @param x       \code{\link{DiffusionMap}} or \code{\link{DPT}}  object
#' @param i,name  Name of a diffusion component 'DCx', 'DPTx', 'Branch' or column from the data
#' @param j       N/A
#' @param ...     ignored
#' 
#' @return The names or data row, see respective generics.
#' 
#' @seealso \link[base]{Extract}, \code{\link[base]{names}} for the generics. \link{DiffusionMap accessors}, \link{DiffusionMap methods}, \link{coercions} for more methods
#' 
#' @examples
#' data(guo)
#' dm <- DiffusionMap(guo)
#' dm$DC1        # A diffusion component
#' dm$Actb       # A gene expression vector
#' dm$num_cells  # Phenotype metadata
#' 
#' dpt <- DPT(dm)
#' dm$Branch
#' dm$DPT1
#' 
#' @name extractions
#' @aliases names.DPT names,DPT-method          $,DPT-method          [[,DPT,character,missing-method
#' names.DiffusionMap names,DiffusionMap-method $,DiffusionMap-method [[,DiffusionMap,character,missing-method
NULL


#' @importFrom methods is
#' @importFrom Biobase featureNames varLabels
#' @name extractions
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
#' @name extractions
#' @export
setMethod('names', 'DPT', function(x) c(paste0('DPT', seq_len(nrow(x))), 'Branch', names(x@dm)))


#' @importFrom methods is
#' @importFrom Biobase exprs featureNames
#' @name extractions
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
#' @name extractions
#' @export
setMethod('[[', c('DPT', 'character', 'missing'), function(x, i, j, ...) {
	if (identical(i, 'dpt')) return(dpt_for_branch(dpt, 1L))  #TODO
	
	num_i <- if (grepl('DPT\\d+', i)) as.integer(sub('DPT(\\d+)', '\\1', i))
	
	if (!is.null(num_i) && 1L <= num_i && num_i <= nrow(x)) {
		x[num_i, ]
	} else if (identical(i, 'Branch') || identical(i, 'branch')) {
		x@branch[, 1L]
	} else x@dm[[i]]
})
#' @name extractions
#' @export
setMethod('[[', c('DPT', 'index', 'index'), function(x, i, j, ...) {
	if (length(i) != 1L || length(j) == 1L)
		stop('Can only extract one element, but i and j were of lengths ', length(i), ' and ', length(j))
	x[i, j, ...]
})


#' @name extractions
#' @export
setMethod('$', 'DiffusionMap', function(x, name) x[[name]])
#' @name extractions
#' @export
setMethod('$', 'DPT', function(x, name) x[[name]])
