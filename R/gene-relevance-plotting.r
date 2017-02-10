#' @include gene-relevance.r
NULL

#' Plot gene relevance or gradient map
#' 
#' \code{plot(gene_relevance, 'Gene')} plots the gradient map of this gene,
#' \code{plot(gene_relevance)} a relevance map of a selection of genes.
#' 
#' @param x            \code{\link{GeneRelevance}} object.
#' @param dm           \code{\link{DiffusionMap}} object.
#' @param y,gene       Gene name or index to create gradient map for. (a number or string)
#' @param iter_smooth  Number of label smoothing iterations to perform on relevance map.
#'                     The higher the more homogenous and the less local structure.
#' @param genes        Genes to based relevance map on or number of genes to use. (vector of strings or one number)
#'                     You can also pass an index into the gene names. (vector of numbers or logicals with length > 1)
#' @param pal          Palette. Either A colormap function or a list of colors
#' @param ...          Passed to \code{plot_gradient_map}/\code{plot_gene_relevance}
#' 
#' @return ggplot2 plot, when plotting a relevance map with a list member $ids containing the IDs used.
#' 
#' @aliases plot.GeneRelevance
#' 
#' @name Gene Relevance plotting
#' @export
setMethod('plot', c('GeneRelevance', 'character'), function(x, y, ...) plot_gradient_map(x@dm, y, ...))
#' @name Gene Relevance plotting
#' @export
setMethod('plot', c('GeneRelevance', 'numeric'),   function(x, y, ...) plot_gradient_map(x@dm, y, ...))

#' @name Gene Relevance plotting
#' @export
setMethod('plot', c('GeneRelevance', 'missing'), function(x, y, ...) plot_gene_relevance(x@dm, 2L, ...))

#' @importFrom ggplot2 ggplot aes geom_point geom_segment scale_colour_gradientn ggtitle
#' @name Gene Relevance plotting
#' @export
plot_gradient_map <- function(dm, gene, ..., pal = cube_helix) {
	if (!is(dm, 'DiffusionMap')) stop('dm needs to be a DiffusionMap object')
	if (length(gene) != 1L) stop('You need to supply a single gene name or index')
	if (is.function(pal)) pal <- pal(12)
	
	relevance_map <- gene_relevance(dm)
	
	all_exprs <- extract_doublematrix(dataset(dm))
	gene_name <- if (is.character(gene)) gene else colnames(all_exprs)[[gene]]
	expr <- all_exprs[, gene]
	partials_norm <- relevance_map@partials_norm[gene, ]
	nn_index <- cbind(seq_along(expr), relevance_map@nn_index)
	
	# Plot gradient vectors
	scatter <- cbind(
		as.data.frame(dm),
		Expression = expr,
		PartialsNorm = partials_norm)
	
	# Select highest vectors in neighbourhoods
	norm_top <- apply(nn_index, 1, function(cell) which.max(partials_norm[cell]) == 1)
	norm_top[sapply(norm_top, length) == 0] <- FALSE
	norm_top <- unlist(norm_top)
	
	dc_var <- .05  # Fraction of overall DC variability
	partials <- lapply(1:2, function(n_dc) {
		dc <- dm[[paste0('DC', n_dc)]][norm_top]
		partials <- relevance_map@partials[gene, norm_top, n_dc]
		# Scale magnitude of partial derivates
		delta <- max(dc, na.rm = TRUE) - min(dc, na.rm = TRUE)
		partials / max(abs(partials), na.rm = TRUE) * dc_var * delta
	})
	
	scatter_top <- cbind(
		scatter[norm_top, ],
		PartialsDC1 = partials[[1]],
		PartialsDC2 = partials[[2]])
	
	ggplot() +
		geom_point(aes(DC1, DC2, colour = Expression), scatter, alpha = 1) + 
		scale_colour_gradientn(colours = pal) + 
		geom_segment(aes(
			x = DC1 - PartialsDC1, xend = DC1 + PartialsDC1,
			y = DC2 - PartialsDC2, yend = DC2 + PartialsDC2,
			alpha = PartialsNorm), scatter_top) +
		ggtitle(gene_name)
}

#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual ggtitle
#' @name Gene Relevance plotting
#' @export
plot_gene_relevance <- function(dm, iter_smooth = 2L, genes = 5L, ..., pal = palette()) {
	if (!is(dm, 'DiffusionMap')) stop('dm needs to be a DiffusionMap object')
	
	relevance_map <- gene_relevance(dm)
	# TODO: non-dm, relevance_map is used below
	partials_norm <- relevance_map@partials_norm
	
	if (is.character(genes)) {
		found <- sapply(genes, function(id) length(grep(id, rownames(partials_norm))) > 0)
		gene_ids <- genes[found]
	} else if (length(genes) == 1L) {
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
	rel_map_data <- cbind(as.data.frame(dm), Gene = as.factor(max_gene))
	
	rel_map <- ggplot(rel_map_data, aes(x = DC1, y = DC2, colour = Gene)) +
		geom_point(alpha = .8) + 
		scale_color_manual(values = pal) +
		ggtitle('Gene relevance map: DC1 vs DC2')
	
	rel_map$ids <- gene_ids
	
	rel_map
}
