#' @include gene-relevance.r
#' @importFrom Biobase featureNames featureNames<-
NULL

#' Gene Relevance methods
#' 
#' \code{featureNames <- ...} Can be used to set the gene names used for plotting
#' (e.g. if the data contains hardly readably gene or transcript IDs)
#' 
#' @seealso \code{\link{gene_relevance}}, \link{Gene Relevance plotting}
#' 
#' @examples
#' data(guo_norm)
#' gr <- gene_relevance(DiffusionMap(guo_norm))
#' featureNames(gr)[[37]] <- 'Id2 (suppresses differentiation)'
#' # now plot it with the changed gene name(s)
#' 
#' @name Gene Relevance methods
#' @export
setMethod('featureNames', 'GeneRelevance', function(object) colnames(object@exprs))

#' @name Gene Relevance methods
#' @export
setMethod('featureNames<-', c('GeneRelevance', 'characterOrFactor'), function(object, value) {
	colnames(object@exprs) <-
		rownames(object@partials_norm) <-
		rownames(object@partials) <-
		value
	object
})
