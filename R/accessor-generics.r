#' destiny generics
#' 
#' destiny provides several generic methods and implements them for the \code{\link{DiffusionMap}} and \code{\link{Sigmas}} classes.
#' 
#' @param object  Object from which to extract or to which to assign a value
#' @param value   Value to assign within an object
#' 
#' @examples
#' data(guo_norm)
#' dm <- DiffusionMap(guo_norm)
#' eigenvalues(dm)
#' eigenvectors(dm)
#' sigmas(dm)
#' optimal_sigma(dm)
#' dataset(dm)
#' distance(dm)
#' 
#' @seealso \link{DiffusionMap methods} and \link{Sigmas} class for implementations
#' 
#' @importFrom methods setGeneric
#' @name destiny generics
#' @rdname destiny-generics
NULL


#' @return \code{eigenvalues} retrieves the numeric eigenvalues
#' @rdname destiny-generics
#' @export
setGeneric('eigenvalues', function(object) standardGeneric('eigenvalues'), valueClass = 'numeric')

#' @rdname destiny-generics
#' @export
setGeneric('eigenvalues<-', function(object, value) standardGeneric('eigenvalues<-'))


#' @return \code{eigenvectors} retrieves the eigenvectors matrix
#' @rdname destiny-generics
#' @export
setGeneric('eigenvectors', function(object) standardGeneric('eigenvectors'), valueClass = 'matrix')

#' @rdname destiny-generics
#' @export
setGeneric('eigenvectors<-', function(object, value) standardGeneric('eigenvectors<-'))


#' @return \code{sigmas} retrieves the \code{\link{Sigmas}} from an object utilizing it as kernel width
#' @rdname destiny-generics
#' @export
setGeneric('sigmas', function(object) standardGeneric('sigmas'), valueClass = 'Sigmas')

#' @rdname destiny-generics
#' @export
setGeneric('sigmas<-', function(object, value) standardGeneric('sigmas<-'))


#' @return \code{dataset} retrieves the data the object was created from
#' @rdname destiny-generics
#' @export
setGeneric('dataset', function(object) standardGeneric('dataset'))

#' @rdname destiny-generics
#' @export
setGeneric('dataset<-', function(object, value) standardGeneric('dataset<-'))


#' @return \code{distance} retrieves the distance metric used to create the object, e.g. \code{euclidean}
#' @rdname destiny-generics
#' @export
setGeneric('distance', function(object) standardGeneric('distance'), valueClass = 'character')

#' @rdname destiny-generics
#' @export
setGeneric('distance<-', function(object, value) standardGeneric('distance<-'))


#' @return \code{optimal_sigma} retrieves the numeric value of the optimal sigma or local sigmas
#' @rdname destiny-generics
#' @export
setGeneric('optimal_sigma', function(object) standardGeneric('optimal_sigma'))

