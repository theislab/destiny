#' @name Gene Relevance plotting
#' @export
setGeneric('plot_gene_relevance', function(coords, exprs, ..., iter_smooth = 2L, n_top = 10L, genes = 5L, dims = 1:2, pal = palette()) {
	standardGeneric('plot_gene_relevance')
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_gene_relevance', c('matrix', 'matrix'), function(coords, exprs, ..., iter_smooth = 2L, n_top = 10L, genes = 5L, dims = 1:2, pal = palette()) {
	plot_gene_relevance_impl(gene_relevance(coords, exprs, dims = seq_len(max(dims))), iter_smooth = iter_smooth, n_top = n_top, genes = genes, dims = dims, pal = pal, ...)
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_gene_relevance', c('DiffusionMap', 'missing'), function(coords, exprs, ..., iter_smooth = 2L, n_top = 10L, genes = 5L, dims = 1:2, pal = palette()) {
	plot_gene_relevance_impl(gene_relevance(coords, dims = seq_len(max(dims))), iter_smooth = iter_smooth, n_top = n_top, genes = genes, dims = dims, pal = pal, ...)
})

#' @name Gene Relevance plotting
#' @export
setMethod('plot_gene_relevance', c('GeneRelevance', 'missing'), function(coords, exprs, ..., iter_smooth = 2L, n_top = 10L, genes = 5L, dims = 1:2, pal = palette()) {
	plot_gene_relevance_impl(coords, iter_smooth = iter_smooth, n_top = n_top, genes = genes, dims = dims, pal = pal, ...)
})

#' @importFrom ggplot2 ggplot aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 ggtitle
#' @importFrom utils head
plot_gene_relevance_impl <- function(relevance_map, ..., iter_smooth, n_top, genes, dims, pal) {
	stopifparams(...)
	relevance_map <- updateObject(relevance_map)
	partials_norm <- relevance_map@partials_norm
	coords <- get_coords(relevance_map, dims)
	if (!is.numeric(iter_smooth) || length(iter_smooth) != 1L) stop('iter_smooth needs to be an integer(1)')
	if (!is.numeric(n_top)       || length(n_top) != 1L      ) stop(      'n_top needs to be an integer(1)')
	
	scores <- NULL
	if (is.character(genes)) {
		found <- sapply(genes, function(id) length(grep(id, colnames(partials_norm))) > 0)
		gene_ids <- genes[found]
	} else if (length(genes) == 1L) {
		n_genes <- min(genes, ncol(relevance_map@exprs), na.rm = TRUE)
		# gene with top n norm for each cell
		genes_max <- if (n_top == 1L) {
			colnames(partials_norm)[apply(partials_norm, 1L, function(cell) which.max(cell))]
		} else {
			genes_ord <- t(apply(partials_norm, 1L, function(cell) {
				cell <- setNames(cell, colnames(partials_norm))
				names(cell[order(cell, decreasing = TRUE)[seq_len(n_top)]])
			}))
			as.vector(genes_ord)
		}
		counts <- as.data.frame(table(genes_max), stringsAsFactors = FALSE)
		counts <- counts[order(counts$Freq, decreasing = TRUE), ]
		n_genes <- min(n_genes, nrow(counts))
		gene_ids <- counts[seq_len(n_genes), 'genes_max']
		scores   <- counts[seq_len(n_genes), 'Freq'] / sum(counts$Freq)
		names(scores) <- gene_ids
	} else {
		gene_ids <- colnames(partials_norm)[genes]
	}
	if (is.function(pal)) pal <- pal(length(gene_ids))
	
	num_top <- min(5L, length(gene_ids))
	top_n_cell_text <- apply(partials_norm, 1L, function(cell) {
		idxs <- head(order(cell, decreasing = TRUE), num_top)
		names <- colnames(partials_norm)[idxs]
		txt <- sprintf('%s. %s (%.3f)', seq_len(num_top), names, cell[idxs])
		paste(txt, collapse = '\n')
	})
	
	# Plot a single map with cells coloured by gene which has 
	# the highest differential norm of all genes considered.
	
	max_gene <-
		if (n_top > 1L) genes_ord[, 1]
		else gene_ids[apply(partials_norm[, gene_ids, drop = FALSE], 1L, function(cell) which.max(cell))]
	# Label smoothing through graph structure
	for (i in seq_len(iter_smooth)) {
		f <-  # for the first smoothing, we can use the full n_top genes
			if (i == 1 && n_top > 1) function(idx) genes_ord[idx, ]
			else function(idx) max_gene[idx]
		max_gene <- apply(relevance_map@nn_index, 1, function(idx_neighbors) {
			max_genes_nn <- unique(f(idx_neighbors))
			max_genes_nn_hist <- sapply(max_genes_nn, function(gene) sum(gene == f(idx_neighbors), na.rm = TRUE))
			print(as.table(max_genes_nn_hist))
			names(max_genes_nn_hist)[which.max(max_genes_nn_hist)]
		})
	}
	# Add more than two DC and return data frame so that user
	# can easily rebuild relevance map on other DC combination than 1 and 2.
	rel_map_data <- cbind(as.data.frame(coords), Gene = factor(max_gene, levels = gene_ids), TopN = top_n_cell_text)
	
	d1 <- colnames(coords)[[1]]
	d2 <- colnames(coords)[[2]]
	rel_map <- ggplot(rel_map_data, aes_string(x = d1, y = d2, colour = 'Gene', text = 'TopN')) +
		geom_point(alpha = .8) +
		geom_rangeframe(colour = par('col')) +
		scale_color_manual(values = pal) +
		ggtitle(sprintf('Gene relevance map')) +
		theme_really_minimal()
	
	rel_map$ids <- gene_ids
	rel_map$scores <- scores
	
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
