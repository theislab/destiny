#' @include diffusionmap.r dpt.r
NULL

#' Coercion methods
#' 
#' Convert a \code{\link{DiffusionMap}} or \code{\link{DPT}} object to other classes
#' 
#' \link[ggplot2]{fortify} is a ggplot2 generic allowing a diffusion map to be used as \code{data} parameter in \link[ggplot2]{ggplot} or \link[ggplot2]{qplot}. 
#' 
#' @param x,model  A \code{\link{DiffusionMap}} or \code{\link{DPT}} object
#' @param row.names  NULL or a character vector giving the row names for the data frame. Missing values are not allowed.
#' @param optional  logical. If TRUE, setting row names and converting column names (to syntactic names: see make.names) is optional.

#' @param ...  Passed to \code{\link[base]{as.data.frame}}
#' @param data  ignored
#' 
#' @return An object of the desired class
#' 
#' @seealso \link{DiffusionMap accessors}, \link{extractions}, \link{DiffusionMap methods} for more methods
#' 
#' @examples
#' library(Biobase)
#' data(guo)
#' dm <- DiffusionMap(guo)
#' classes <- vapply(as.data.frame(dm), class, character(1L))
#' stopifnot(all(classes[paste0('DC', 1:20)] == 'numeric'))
#' stopifnot(all(classes[featureNames(guo) ] == 'numeric'))
#' stopifnot(all(classes[   varLabels(guo) ] == c('factor', 'integer')))
#' 
#' @aliases
#' as.data.frame.DiffusionMap as.data.frame,DiffusionMap-method fortify.DiffusionMap
#' as.data.frame.DPT          as.data.frame,DPT-method          fortify.DPT
#'     as.matrix.DPT              as.matrix,DPT-method
#' 
#' @importFrom BiocGenerics as.data.frame
#' @name coercions
#' @include diffusionmap.r
NULL


#' @name coercions
#' @export
setMethod('as.data.frame', 'DiffusionMap', function(x, row.names = NULL, optional = FALSE, ...) {
	cbind(as.data.frame(eigenvectors(x), row.names, optional, ...),
				as.data.frame(dataset(x),      row.names, optional, ...))
})


#' @usage fortify.DiffusionMap(model, data, ...)
#' 
#' @importFrom BiocGenerics as.data.frame
#' @importFrom Biobase as.data.frame.ExpressionSet
#' @name coercions
#' @export fortify.DiffusionMap
fortify.DiffusionMap <- function(model, data, ...) as.data.frame(model, ...)


#' @name coercions
#' @export
setMethod('as.data.frame', 'DPT', function(x, row.names = NULL, optional = FALSE, ...) cbind(
	as.data.frame(cbind(x[], x@branch, x@tips), row.names = row.names, optional = optional, ...),
	as.data.frame(x@dm,                         row.names = row.names, optional = optional, ...)))


#' @usage fortify.DPT(model, data, ...)
#' 
#' @importFrom BiocGenerics as.data.frame
#' @importFrom Biobase as.data.frame.ExpressionSet
#' @name coercions
#' @export fortify.DPT
fortify.DPT <- function(model, data, ...) as.data.frame(model, ...)

#' @name coercions
#' @export
setMethod('as.matrix', 'DPT', function(x, ...) x[])
