#' DiffusionMap coercion methods
#' 
#' Convert a \code{\link{DiffusionMap}} to other classes
#' 
#' \link[ggplot2]{fortify} is a ggplot2 generic allowing a diffusion map to be used as \code{data} parameter in \link[ggplot2]{ggplot} or \link[ggplot2]{qplot}. 
#' 
#' @param x,model  A \code{\link{DiffusionMap}}
#' @param data  ignored
#' @param row.names  NULL or a character vector giving the row names for the data frame. Missing values are not allowed.
#' @param optional  logical. If TRUE, setting row names and converting column names (to syntactic names: see make.names) is optional.
#' @param ...  Other parameters
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
#' @aliases as.data.frame.DiffusionMap as.data.frame,DiffusionMap-method fortify.DiffusionMap
#' @importFrom BiocGenerics as.data.frame
#' @name DiffusionMap coercions
#' @include diffusionmap.r
NULL


#' @name DiffusionMap coercions
#' @export
setMethod('as.data.frame', 'DiffusionMap', function(x, row.names = NULL, optional = FALSE, ...) {
	cbind(as.data.frame(eigenvectors(x), row.names, optional, ...),
				as.data.frame(dataset(x),      row.names, optional, ...))
})


#' @name DiffusionMap coercions
#' @usage fortify.DiffusionMap(model, data, ...)
#' @export fortify.DiffusionMap
fortify.DiffusionMap <- function(model, data, ...) BiocGenerics::as.data.frame(model)
