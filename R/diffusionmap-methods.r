#' DiffusionMap methods
#' 
#' Methods for external operations on diffusion maps
#' 
#' @param x,object  A \code{\link{DiffusionMap}}
#' 
#' @return The \code{DiffusionMap} object (\code{print}), or NULL (\code{show}), invisibly
#' 
#' @seealso \link{DiffusionMap accessors}, \link{extractions}, \link{coercions} for more methods
#' 
#' @examples
#' data(guo)
#' dm <- DiffusionMap(guo)
#' print(dm)
#' show(dm)
#' 
#' @aliases DiffusionMap-methods
#' print.DiffusionMap        show.DiffusionMap       
#' print,DiffusionMap-method show,DiffusionMap-method
#' @name DiffusionMap methods
NULL

#' @importFrom utils str
#' 
#' @name DiffusionMap methods
#' @export
setMethod('print', 'DiffusionMap', function(x) {
	cat(sprintf('DiffusionMap (%s Diffusion components and %s samples)\n', length(eigenvalues(x)), nrow(eigenvectors(x))))
	cat('eigenvalues:   '); str(eigenvalues(x))
	cat('eigenvectors:  '); str(structure(eigenvectors(x), dimnames = NULL))
	cat('  ..colnames:  '); str(colnames(eigenvectors(x)), vec.len = 4)
	cat('optimal_sigma: '); str(optimal_sigma(x))
	cat('distance:      '); str(distance(x))
	invisible(x)
})

#' @name DiffusionMap methods
#' @export
setMethod('show', 'DiffusionMap', function(object) {
	print(object)
	invisible()
})
