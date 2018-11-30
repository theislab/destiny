#' @include gene-relevance.r
NULL

#' Plot gene relevance or differential map
#' 
#' \code{plot(gene_relevance, 'Gene')} plots the differential map of this/these gene(s),
#' \code{plot(gene_relevance)} a relevance map of a selection of genes.
#' Alternatively, you can use \code{plot_differential_map} or \code{plot_gene_relevance} on a \code{\link{GeneRelevance}} or \code{\link{DiffusionMap}} object, or with two matrices.
#' 
#' @param x            \code{\link{GeneRelevance}} object.
#' @param y,gene       Gene name(s) or index/indices to create differential map for. (integer or character)
#' @param coords       A \code{\link{DiffusionMap}}/\code{\link{GeneRelevance}} object or a cells \eqn{\times} dims \code{\link{matrix}}.
#' @param exprs        An cells \eqn{\times} genes \code{\link{matrix}}. Only provide if \code{coords} is a matrix.
#' @param ...          Passed to \code{plot_differential_map}/\code{plot_gene_relevance}.
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
#' plot(gr, 'Fgf4')  # or plot_differential_map(dm, 'Fgf4')
#' 
#' guo_norm_mat <- t(Biobase::exprs(guo_norm))
#' pca <- prcomp(guo_norm_mat)$x
#' plot_gene_relevance(pca, guo_norm_mat, dims = 2:3)
#' plot_differential_map(pca, guo_norm_mat, gene = c('Fgf4', 'Nanog'))
#' 
#' @aliases
#'   plot.GeneRelevance
#'   plot,GeneRelevance,character-method
#'   plot,GeneRelevance,numeric-method
#'   plot,GeneRelevance,missing-method
#'   plot_differential_map
#'   plot_differential_map,matrix,matrix-method
#'   plot_differential_map,DiffusionMap,missing-method
#'   plot_differential_map,GeneRelevance,missing-method
#'   plot_gene_relevance
#'   plot_gene_relevance,matrix,matrix-method
#'   plot_gene_relevance,DiffusionMap,missing-method
#'   plot_gene_relevance,GeneRelevance,missing-method
#' 
#' @name Gene Relevance plotting
#' @export
setMethod('plot', c('GeneRelevance', 'character'), function(x, y, ...) plot_differential_map(x, gene = y, ...))
#' @name Gene Relevance plotting
#' @export
setMethod('plot', c('GeneRelevance', 'numeric'),   function(x, y, ...) plot_differential_map(x, gene = y, ...))

#' @name Gene Relevance plotting
#' @export
setMethod('plot', c('GeneRelevance', 'missing'), function(x, y, ...) plot_gene_relevance(x, ...))


# plot_differential_map -------------------------------------------------------------------------------------------------------


#' @name Gene Relevance plotting
#' @export
setGeneric('plot_differential_map', function(coords, exprs, ..., gene, dims = 1:2, pal = cube_helix, faceter = facet_wrap(~ Gene)) {
	standardGeneric('plot_differential_map')
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_differential_map', c('matrix', 'matrix'), function(coords, exprs, ..., gene, dims = 1:2, pal = cube_helix, faceter = facet_wrap(~ Gene)) {
	plot_differential_map_impl(gene_relevance(coords, exprs, dims = seq_len(max(dims))), genes = gene, dims = dims, pal = pal, faceter = faceter)
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_differential_map', c('DiffusionMap', 'missing'), function(coords, exprs, ..., gene, dims = 1:2, pal = cube_helix, faceter = facet_wrap(~ Gene)) {
	plot_differential_map_impl(gene_relevance(coords, dims = seq_len(max(dims))), genes = gene, dims = dims, pal = pal, faceter = faceter)
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_differential_map', c('GeneRelevance', 'missing'), function(coords, exprs, ..., gene, dims = 1:2, pal = cube_helix, faceter = facet_wrap(~ Gene)) {
	plot_differential_map_impl(coords, genes = gene, dims = dims, pal = pal, faceter = faceter)
})

differential_map <- function(relevance_map, genes = NULL, dims = 1:2, all = FALSE) {
	relevance_map <- updateObject(relevance_map)
	
	all_dims <- get_dim_range(relevance_map@partials, 3L, dims)
	if (!all(dims %in% all_dims)) stop(
		'The relevance map contains only the dimensions ', paste(all_dims, collapse = ', '),
		', not ', paste(setdiff(dims, all_dims), collapse = ', '))
	
	genes_existing <- colnames(relevance_map@partials_norm)
	if (is.null(genes)) genes <- genes_existing
	else {
		genes_missing <- is.na(match(genes, genes_existing))
		if (any(genes_missing)) stop(
			'The dataset used for the relevance map does not contain gene(s) ', paste(genes[genes_missing], collapse = ', '),
			'. Did you mean ', paste(agrep(genes[genes_missing], genes_existing, value = TRUE), collapse = ', '), '?')
	}
	
	exprs <- relevance_map@exprs
	coords <- get_coords(relevance_map, dims)
	
	partials_norms <- relevance_map@partials_norm[, genes, drop = FALSE]
	nn_index <- cbind(seq_len(nrow(exprs)), relevance_map@nn_index)
	
	# Plot differential vectors
	scatters <- do.call(rbind, lapply(genes, function(g) {
		cbind(
			as.data.frame(coords),
			Expression = exprs[, g],
			PartialsNorm = partials_norms[, g],
			Gene = g)
	}))
	
	scatters_top <- do.call(rbind, lapply(genes, function(g) {
		norm_top <- if (all) seq_len(nrow(exprs)) else {
			# Select highest vectors in neighbourhoods
			norm_top <- apply(nn_index, 1, function(cell) which.max(partials_norms[cell, g]) == 1)
			norm_top[sapply(norm_top, length) == 0] <- FALSE
			norm_top <- unlist(norm_top)
		}
		
		d_var <- .05  # Fraction of overall dimension variability
		partials <- lapply(seq_len(length(dims)), function(d) {
			dc <- coords[norm_top, d]
			partials <- relevance_map@partials[norm_top, g, d]
			# Scale magnitude of partial derivates
			delta <- diff(rev(range(dc, na.rm = TRUE)))
			partials / max(abs(partials), na.rm = TRUE) * d_var * delta
		})
		
		scatter <- subset(scatters, scatters$Gene == g)
		cbind(
			scatter[norm_top, ],
			Angle     = atan(partials[[2]]   / partials[[1]]  ),
			Magnitude = sqrt(partials[[1]]^2 + partials[[2]]^2)
		)
	}))
	
	if (all) scatters_top
	else list(scatters = scatters, scatters_top = scatters_top)
}

#' @importFrom ggplot2 ggplot aes aes_string
#' @importFrom ggplot2 geom_point geom_spoke
#' @importFrom ggplot2 scale_colour_gradientn
#' @importFrom ggplot2 ggtitle facet_wrap
#' @importFrom grid arrow unit
plot_differential_map_impl <- function(relevance_map, ..., genes, dims, pal, faceter) {
	if (is.function(pal)) pal <- pal(12)
	dtm <- differential_map(relevance_map, genes, dims)
	coords <- get_coords(relevance_map, dims)
	gene_names <- if (is.character(genes)) genes else colnames(relevance_map@exprs)[genes]
	
	d1 <- colnames(coords)[[1]]
	d2 <- colnames(coords)[[2]]
	gg <- ggplot() +
		geom_point(aes_string(d1, d2, colour = 'Expression'), dtm$scatters, alpha = 1) + 
		scale_colour_gradientn(colours = pal) + 
		geom_spoke(
			aes_string(
				x = 'DC1', y = 'DC2',
				angle = 'Angle', radius = 'Magnitude',
				alpha = 'PartialsNorm'),
			dtm$scatters_top,
			arrow = arrow(length = unit(.01, 'npc')))
	
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

#' @importFrom ggplot2 ggplot aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 ggtitle
#' @importFrom utils head
plot_gene_relevance_impl <- function(relevance_map, ..., iter_smooth, genes, dims, pal) {
	relevance_map <- updateObject(relevance_map)
	partials_norm <- relevance_map@partials_norm
	coords <- get_coords(relevance_map, dims)
	
	if (is.character(genes)) {
		found <- sapply(genes, function(id) length(grep(id, colnames(partials_norm))) > 0)
		gene_ids <- genes[found]
	} else if (length(genes) == 1L) {
		n_genes <- min(genes, ncol(relevance_map@exprs), na.rm = TRUE)
		# gene with max norm for each cell
		genes_max <- colnames(partials_norm)[apply(partials_norm, 1L, function(cell) which.max(cell))]
		counts <- as.data.frame(table(genes_max), stringsAsFactors = FALSE)
		n_genes <- min(n_genes, nrow(counts))
		gene_ids <- counts[order(counts$Freq, decreasing = TRUE)[1:n_genes], 'genes_max']
	} else {
		gene_ids <- colnames(partials_norm)[genes]
	}
	if (is.function(pal)) pal <- pal(length(gene_ids))
	
	num_top <- min(5L, length(gene_ids))
	top_n <- apply(partials_norm, 1L, function(cell) {
		idxs <- head(order(cell, decreasing = TRUE), num_top)
		names <- colnames(partials_norm)[idxs]
		txt <- sprintf('%s. %s (%.3f)', seq_len(num_top), names, cell[idxs])
		paste(txt, collapse = '\n')
	})
	
	# Plot a single map with cells coloured by gene which has 
	# the highest differential norm of all genes considered.
	
	max_gene_idx <- apply(partials_norm[, gene_ids, drop = FALSE], 1L, which.max)
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
	rel_map_data <- cbind(as.data.frame(coords), Gene = factor(max_gene, levels = gene_ids), TopN = top_n)
	
	d1 <- colnames(coords)[[1]]
	d2 <- colnames(coords)[[2]]
	rel_map <- ggplot(rel_map_data, aes_string(x = d1, y = d2, colour = 'Gene', text = 'TopN')) +
		geom_point(alpha = .8) + 
		scale_color_manual(values = pal) +
		ggtitle(sprintf('Gene relevance map'))
	
	rel_map$ids <- gene_ids
	
	rel_map
}

get_coords <- function(relevance_map, dims) {
	coords <- relevance_map@coords[, dims, drop = FALSE]
	if (is.null(colnames(coords)))
		colnames(coords) <- paste0('Dim', dims)
	coords
}

get_dim_range <- function(arr, d, dims) {
	if (is.character(dims)) dimnames(arr)[[d]]
	else seq_len(dim(arr)[[d]])
}
