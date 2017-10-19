#' @include gene-relevance.r
NULL

#' Plot gene relevance or gradient map
#' 
#' \code{plot(gene_relevance, 'Gene')} plots the gradient map of this/these gene(s),
#' \code{plot(gene_relevance)} a relevance map of a selection of genes.
#' Alternatively, you can use \code{plot_gradient_map} or \code{plot_gene_relevance} on a \code{\link{GeneRelevance}} or \code{\link{DiffusionMap}} object, or with two matrices.
#' 
#' @param x            \code{\link{GeneRelevance}} object.
#' @param y,gene       Gene name(s) or index/indices to create gradient map for. (integer or character)
#' @param coords       A \code{\link{DiffusionMap}}/\code{\link{GeneRelevance}} object or a cells \eqn{\times} dims \code{\link{matrix}}.
#' @param exprs        An cells \eqn{\times} genes \code{\link{matrix}}. Only provide if \code{coords} is a matrix.
#' @param ...          Passed to \code{plot_gradient_map}/\code{plot_gene_relevance}.
#' @param iter_smooth  Number of label smoothing iterations to perform on relevance map.
#'                     The higher the more homogenous and the less local structure.
#' @param genes        Genes to based relevance map on or number of genes to use. (vector of strings or one number)
#'                     You can also pass an index into the gene names. (vector of numbers or logicals with length > 1)
#' @param dims         Names or indices of dimensions to plot. When not plotting a \code{\link{GeneRelevance}} object, the relevance for the dimensions \code{1:max(dims)} will be calculated.
#' @param pal          Palette. Either A colormap function or a list of colors.
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
#' plot(gr, 'Fgf4')  # or plot_gradient_map(dm, 'Fgf4')
#' 
#' guo_norm_mat <- t(Biobase::exprs(guo_norm))
#' pca <- prcomp(guo_norm_mat)$x
#' plot_gene_relevance(pca, guo_norm_mat, dims = 2:3)
#' plot_gradient_map(pca, guo_norm_mat, gene = c('Fgf4', 'Nanog'))
#' 
#' @aliases
#'   plot.GeneRelevance
#'   plot,GeneRelevance,character-method
#'   plot,GeneRelevance,numeric-method
#'   plot,GeneRelevance,missing-method
#'   plot_gradient_map
#'   plot_gradient_map,matrix,matrix-method
#'   plot_gradient_map,DiffusionMap,missing-method
#'   plot_gradient_map,GeneRelevance,missing-method
#'   plot_gene_relevance
#'   plot_gene_relevance,matrix,matrix-method
#'   plot_gene_relevance,DiffusionMap,missing-method
#'   plot_gene_relevance,GeneRelevance,missing-method
#' 
#' @name Gene Relevance plotting
#' @export
setMethod('plot', c('GeneRelevance', 'character'), function(x, y, ...) plot_gradient_map(x, gene = y, ...))
#' @name Gene Relevance plotting
#' @export
setMethod('plot', c('GeneRelevance', 'numeric'),   function(x, y, ...) plot_gradient_map(x, gene = y, ...))

#' @name Gene Relevance plotting
#' @export
setMethod('plot', c('GeneRelevance', 'missing'), function(x, y, ...) plot_gene_relevance(x, ...))


# plot_gradient_map -------------------------------------------------------------------------------------------------------


#' @name Gene Relevance plotting
#' @export
setGeneric('plot_gradient_map', function(coords, exprs, ..., gene, dims = 1:2, pal = cube_helix, faceter = facet_wrap(~ Gene)) standardGeneric('plot_gradient_map'))

#' @name Gene Relevance plotting
#' @export
setMethod('plot_gradient_map', c('matrix', 'matrix'), function(coords, exprs, ..., gene, dims = 1:2, pal = cube_helix, faceter = facet_wrap(~ Gene)) {
	plot_gradient_map_impl(gene_relevance(coords, exprs, dims = seq_len(max(dims))), genes = gene, dims = dims, pal = pal, faceter = faceter)
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_gradient_map', c('DiffusionMap', 'missing'), function(coords, exprs, ..., gene, dims = 1:2, pal = cube_helix, faceter = facet_wrap(~ Gene)) {
	plot_gradient_map_impl(gene_relevance(coords, dims = seq_len(max(dims))), genes = gene, dims = dims, pal = pal, faceter = faceter)
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_gradient_map', c('GeneRelevance', 'missing'), function(coords, exprs, ..., gene, dims = 1:2, pal = cube_helix, faceter = facet_wrap(~ Gene)) {
	plot_gradient_map_impl(coords, genes = gene, dims = dims, pal = pal, faceter = faceter)
})

#' @importFrom ggplot2 ggplot aes aes_string geom_point geom_segment scale_colour_gradientn ggtitle facet_wrap
plot_gradient_map_impl <- function(relevance_map, ..., genes, dims, pal, faceter) {
	if (missing(genes)) stop('You need to supply gene name(s) or index/indices')
	if (is.function(pal)) pal <- pal(12)
	
	exprs <- relevance_map@exprs
	coords <- get_coords(relevance_map, dims)
	
	gene_names <- if (is.character(genes)) genes else colnames(exprs)[genes]
	partials_norms <- relevance_map@partials_norm[genes, , drop = FALSE]
	nn_index <- cbind(seq_len(nrow(exprs)), relevance_map@nn_index)
	
	# Plot gradient vectors
	scatters <- do.call(rbind, lapply(genes, function(g) {
		cbind(
			as.data.frame(coords),
			Expression = exprs[, g],
			PartialsNorm = partials_norms[g, ],
			Gene = g)
	}))
	
	scatters_top <- do.call(rbind, lapply(genes, function(g) {
		# Select highest vectors in neighbourhoods
		norm_top <- apply(nn_index, 1, function(cell) which.max(partials_norms[g, cell]) == 1)
		norm_top[sapply(norm_top, length) == 0] <- FALSE
		norm_top <- unlist(norm_top)
		
		d_var <- .05  # Fraction of overall dimension variability
		partials <- lapply(seq_len(length(dims)), function(d) {
			dc <- coords[norm_top, d]
			partials <- relevance_map@partials[g, norm_top, d]
			# Scale magnitude of partial derivates
			delta <- diff(rev(range(dc, na.rm = TRUE)))
			partials / max(abs(partials), na.rm = TRUE) * d_var * delta
		})
		
		scatter <- subset(scatters, Gene == g)
		D1 <- scatter[norm_top, 1]
		D2 <- scatter[norm_top, 2]
		cbind(
			scatter[norm_top, ],
			D1start = D1 - partials[[1]], D1end = D1 + partials[[1]],
			D2start = D2 - partials[[2]], D2end = D2 + partials[[2]],
			Gene = g)
	}))
	
	d1 <- colnames(coords)[[1]]
	d2 <- colnames(coords)[[2]]
	gg <- ggplot() +
		geom_point(aes_string(d1, d2, colour = 'Expression'), scatters, alpha = 1) + 
		scale_colour_gradientn(colours = pal) + 
		geom_segment(aes_string(
			x = 'D1start', xend = 'D1end',
			y = 'D2start', yend = 'D2end',
			alpha = 'PartialsNorm'), scatters_top)
	
	if (length(genes) > 1) gg + faceter else gg + ggtitle(gene_names)
}


# plot_gene_relevance -----------------------------------------------------------------------------------------------------

#' @name Gene Relevance plotting
#' @export
setGeneric('plot_gene_relevance', function(coords, exprs, ..., iter_smooth = 2L, genes = 5L, dims = 1:2, pal = palette()) standardGeneric('plot_gene_relevance'))

#' @name Gene Relevance plotting
#' @export
setMethod('plot_gene_relevance', c('matrix', 'matrix'), function(coords, exprs, ..., iter_smooth = 2L, genes = 5L, dims = 1:2, pal = palette()) {
	plot_gene_relevance_impl(gene_relevance(coords, exprs, dims = seq_len(max(dims))), iter_smooth = iter_smooth, genes = genes, dims = dims, pal = pal)
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_gene_relevance', c('DiffusionMap', 'missing'), function(coords, exprs, ..., iter_smooth = 2L, genes = 5L, dims = 1:2, pal = palette()) {
	plot_gene_relevance_impl(gene_relevance(coords, dims = seq_len(max(dims))), iter_smooth = iter_smooth, genes = genes, dims = dims, pal = pal)
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_gene_relevance', c('GeneRelevance', 'missing'), function(coords, exprs, ..., iter_smooth = 2L, genes = 5L, dims = 1:2, pal = palette()) {
	plot_gene_relevance_impl(coords, iter_smooth = iter_smooth, genes = genes, dims = dims, pal = pal)
})

#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_manual ggtitle
plot_gene_relevance_impl <- function(relevance_map, ..., iter_smooth, genes, dims, pal) {
	partials_norm <- relevance_map@partials_norm
	
	coords <- get_coords(relevance_map, dims)
	
	if (is.character(genes)) {
		found <- sapply(genes, function(id) length(grep(id, rownames(partials_norm))) > 0)
		gene_ids <- genes[found]
	} else if (length(genes) == 1L) {
		genes <- min(genes, ncol(relevance_map@exprs))
		# Select most often occurring genes with maximal norm of gradients at a cell.
		max_gene_ids <- rownames(partials_norm)[unlist(apply(partials_norm, 2, function(cell) which.max(cell)))]
		max_gene_ids_hist <- sapply(unique(max_gene_ids), function(id) sum(max_gene_ids == id, na.rm = TRUE) )
		gene_ids <- names(max_gene_ids_hist)[sort(max_gene_ids_hist, decreasing = TRUE, index.return = TRUE)$ix[1:genes]]
	} else {
		gene_ids <- rownames(partials_norm)[genes]
	}
	if (is.function(pal)) pal <- pal(length(gene_ids))
	
	# Plot a single map with cells coloured by gene which has 
	# the highest gradient norm of all genes considered.
	
	max_gene_idx <- apply(partials_norm[gene_ids, , drop = FALSE], 2, which.max)
	max_gene_idx[sapply(max_gene_idx, length) == 0] <- NA
	max_gene <- gene_ids[unlist(max_gene_idx)]
	# Label smoothing through graph structure
	for (i in seq_len(iter_smooth)) {
		max_gene <- apply(relevance_map@nn_index, 1, function(cell) {
			max_genes_nn <- unique(max_gene[cell])
			max_genes_nn_hist <- sapply(max_genes_nn, function(gene) sum(gene == max_gene[cell], na.rm = TRUE))
			names(max_genes_nn_hist)[which.max(max_genes_nn_hist)]
		})
	}
	# Add more than two DC and return data frame so that user
	# can easily rebuild relevance map on other DC combination than 1 and 2.
	rel_map_data <- cbind(as.data.frame(coords), Gene = as.factor(max_gene))
	
	d1 <- colnames(coords)[[1]]
	d2 <- colnames(coords)[[2]]
	rel_map <- ggplot(rel_map_data, aes_string(x = d1, y = d2, colour = 'Gene')) +
		geom_point(alpha = .8) + 
		scale_color_manual(values = pal) +
		ggtitle(sprintf('Gene relevance map'))
	
	rel_map$ids <- gene_ids
	
	rel_map
}

get_coords <- function(relevance_map, dims) {
	coords <- relevance_map@coords[, dims, drop = FALSE]
	if (is.null(colnames(coords)))
		colnames(coords) <- paste('Dimension', dims)
	coords
}
