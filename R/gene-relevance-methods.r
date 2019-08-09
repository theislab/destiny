#' @include gene-relevance.r
#' @importFrom Biobase featureNames featureNames<-
NULL

#' Gene Relevance methods
#' 
#' \code{featureNames <- ...} can be used to set the gene names used for plotting
#' (e.g. if the data contains hardly readably gene or transcript IDs).
#' \code{dataset} gets the expressions used for the gene relevance calculations,
#' and \code{distance} the distance measure.
#' 
#' @param x,object  \code{GeneRelevance} object
#' @param value     A text vector (\code{\link{character}} or \code{\link{factor}})
#' 
#' @return
#' \code{dataset}, \code{distance}, and \code{featureNames} return the stored properties.
#' The other methods return a \code{GeneRelevance} object (\code{print}, \code{... <- ...}),
#' or NULL (\code{show}), invisibly
#' 
#' @seealso \code{\link{gene_relevance}}, \link{Gene Relevance plotting}
#' 
#' @examples
#' data(guo_norm)
#' dm <- DiffusionMap(guo_norm)
#' gr <- gene_relevance(dm)
#' stopifnot(distance(gr) == distance(dm))
#' featureNames(gr)[[37]] <- 'Id2 (suppresses differentiation)'
#' # now plot it with the changed gene name(s)
#' 
#' @aliases featureNames.GeneRelevance dataset.GeneRelevance
#' @name Gene Relevance methods
#' @rdname Gene-Relevance-methods
NULL


#' @importFrom methods is
#' @importFrom utils str
#' 
#' @rdname Gene-Relevance-methods
#' @export
setMethod('print', 'GeneRelevance', function(x) {
	d <- dataset(x)
	cat(sprintf('GeneRelevance (%s genes, %s reduced dimensions, and %s observations)\n', ncol(x@exprs), ncol(x@coords), nrow(x@exprs)))
	cat('is:      ')
	if (is(d, 'Matrix')) cat(sprintf('%s%s%s %s (%s)\n', nrow(d), sym_times, ncol(d), class(d)[[1L]], mode(d@x)))
	else str(structure(d, dimnames = NULL))
	cat('featureNames: '); str(featureNames(x))
	invisible(x)
})

#' @importFrom methods show
#' 
#' @rdname Gene-Relevance-methods
#' @export
setMethod('show', 'GeneRelevance', function(object) {
	print(object)
	invisible()
})


#' @rdname Gene-Relevance-methods
#' @export
setMethod('featureNames', 'GeneRelevance', function(object) colnames(object@exprs))

#' @rdname Gene-Relevance-methods
#' @export
setMethod('featureNames<-', c('GeneRelevance', 'characterOrFactor'), function(object, value) {
	object <- updateObject(object)
	colnames(object@exprs) <-
		colnames(object@partials_norm) <-
		colnames(object@partials) <-
		value
	object
})


#' @rdname Gene-Relevance-methods
#' @export
setMethod('dataset', 'GeneRelevance', function(object) object@exprs)

#' @rdname Gene-Relevance-methods
#' @export
setMethod('dataset<-', 'GeneRelevance', function(object, value) {
	object <- updateObject(object)
	object@exprs <- value
	validObject(object)
	object
})


#' @rdname Gene-Relevance-methods
#' @export
setMethod('distance', 'GeneRelevance', function(object) object@distance)

#' @rdname Gene-Relevance-methods
#' @export
setMethod('distance<-', 'GeneRelevance', function(object, value) {
	object@distance <- value
	validObject(object)
	object
})
