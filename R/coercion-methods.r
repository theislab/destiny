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
#' @seealso \link{DiffusionMap accessors}, \link{DiffusionMap extraction}, \link{DiffusionMap methods} for more methods
#' 
#' @examples
#' data(guo)
#' dm <- DiffusionMap(guo)
#' dm$DC1   # A diffusion component
#' dm$Actb  # A gene expression
#' 
#' @aliases
#' as.data.frame.DiffusionMap as.data.frame,DiffusionMap-method fortify.DiffusionMap
#' as.data.frame.DPT          as.data.frame,DPT-method          fortify.DPT
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
	data.frame(
		Branch = x@branch,
		Parent = x@parent,
		Tip    = lWhich(x@tips, len = length(x@branch)),
		row.names = row.names),
	as.data.frame(x@dpt, row.names = row.names, optional = optional, ...),
	as.data.frame(x@dm,  row.names = row.names, optional = optional, ...)))


#' @usage fortify.DPT(model, data, ...)
#' 
#' @importFrom BiocGenerics as.data.frame
#' @importFrom Biobase as.data.frame.ExpressionSet
#' @name coercions
#' @export fortify.DPT
fortify.DPT <- function(model, data, ...) as.data.frame(model, ...)
