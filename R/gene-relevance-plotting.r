#' @include gene-relevance.r
NULL

#' Plot gene relevance or differential map
#' 
#' \code{plot(gene_relevance, 'Gene')} plots the differential map of this/these gene(s),
#' \code{plot(gene_relevance)} a relevance map of a selection of genes.
#' Alternatively, you can use \code{plot_differential_map} or \code{plot_gene_relevance} on a \code{\link{GeneRelevance}} or \code{\link{DiffusionMap}} object, or with two matrices.
#' 
#' @param x            \code{\link{GeneRelevance}} object.
#' @param y            Gene name(s) or index/indices to create differential map for. (integer or character)
#' @param coords       A \code{\link{DiffusionMap}}/\code{\link{GeneRelevance}} object or a cells \eqn{\times} dims \code{\link{matrix}}.
#' @param exprs        An cells \eqn{\times} genes \code{\link{matrix}}. Only provide if \code{coords} is a matrix.
#' @param ...          Passed to \code{plot_differential_map}/\code{plot_gene_relevance}.
#' @param iter_smooth  Number of label smoothing iterations to perform on relevance map.
#'                     The higher the more homogenous and the less local structure.
#' @param n_top        Number the top n genes per cell count towards the score defining which genes to return and plot in the relevance map.
#' @param genes        Genes to base relevance map on (vector of strings).
#'                     You can also pass an index into the gene names (vector of numbers or logicals with length > 1).
#'                     The default NULL means all genes.
#' @param dims         Names or indices of dimensions to plot. When not plotting a \code{\link{GeneRelevance}} object, the relevance for the dimensions \code{1:max(dims)} will be calculated.
#' @param pal          Palette. Either A colormap function or a list of colors.
#' @param col_na       Color for cells that end up with no most relevant gene.
#' @param limit        Limit the amount of displayed gene labels to the amount of available colors in \code{pal}?
#' @param faceter      A ggplot faceter like \code{\link[ggplot2]{facet_wrap}(~ Gene)}.
#' 
#' @return ggplot2 plot, when plotting a relevance map with a list member \code{$ids} containing the gene IDs used.
#' 
#' @seealso \code{\link{gene_relevance}}, \link{Gene Relevance methods}
#' 
#' @examples
#' data(guo_norm)
#' dm <- DiffusionMap(guo_norm)
#' gr <- gene_relevance(dm)
#' plot(gr)          # or plot_gene_relevance(dm)
#' plot(gr, 'Fgf4')  # or plot_differential_map(dm, 'Fgf4')
#' 
#' guo_norm_mat <- t(Biobase::exprs(guo_norm))
#' pca <- prcomp(guo_norm_mat)$x
#' plot_gene_relevance(pca, guo_norm_mat, dims = 2:3)
#' plot_differential_map(pca, guo_norm_mat, genes = c('Fgf4', 'Nanog'))
#' 
#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot', c('GeneRelevance', 'character'), function(x, y, ...) plot_differential_map(x, genes = y, ...))
#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot', c('GeneRelevance', 'numeric'),   function(x, y, ...) plot_differential_map(x, genes = y, ...))

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot', c('GeneRelevance', 'missing'), function(x, y, ...) plot_gene_relevance(x, ...))
