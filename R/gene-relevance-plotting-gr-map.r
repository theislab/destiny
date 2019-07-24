#' @rdname Gene-Relevance-plotting
#' @export
setGeneric('plot_gene_relevance', function(coords, exprs, ..., iter_smooth = 2L, n_top = 10L, genes = NULL, dims = 1:2, pal = palette(), col_na = 'grey', limit = TRUE) {
	standardGeneric('plot_gene_relevance')
})

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot_gene_relevance', c('matrix', 'matrix'), function(coords, exprs, ..., iter_smooth, n_top, genes, dims, pal, col_na, limit) {
	plot_gene_relevance_impl(gene_relevance(coords, exprs, dims = seq_len(max(dims))), iter_smooth = iter_smooth, n_top = n_top, genes = genes, dims = dims, pal = pal, col_na = col_na, limit = limit, ...)
})

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot_gene_relevance', c('DiffusionMap', 'missing'), function(coords, exprs, ..., iter_smooth, n_top, genes, dims, pal, col_na, limit) {
	plot_gene_relevance_impl(gene_relevance(coords, dims = seq_len(max(dims))), iter_smooth = iter_smooth, n_top = n_top, genes = genes, dims = dims, pal = pal, col_na = col_na, limit = limit, ...)
})

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot_gene_relevance', c('GeneRelevance', 'missing'), function(coords, exprs, ..., iter_smooth, n_top, genes, dims, pal, col_na, limit) {
	plot_gene_relevance_impl(coords, iter_smooth = iter_smooth, n_top = n_top, genes = genes, dims = dims, pal = pal, col_na = col_na, limit = limit, ...)
})

#' @importFrom ggplot2 ggplot aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 ggtitle
#' @importFrom Biobase featureNames
#' @importFrom utils head
plot_gene_relevance_impl <- function(relevance_map, ..., iter_smooth, n_top, genes, dims, pal, col_na, limit) {
	stopifparams(...)
	relevance_map <- updateObject(relevance_map)
	partials_norm <- relevance_map@partials_norm
	coords <- get_coords(relevance_map, dims)
	if (!is.numeric(iter_smooth) || length(iter_smooth) != 1L) stop('iter_smooth needs to be an integer(1)')
	if (!is.numeric(n_top)       || length(n_top) != 1L      ) stop(      'n_top needs to be an integer(1)')
	
	all_genes <- featureNames(relevance_map)
	if (is.null(genes)) {
		genes <- all_genes
	} else if (is.character(genes)) {
		found <- sapply(genes, function(id) length(grep(id, all_genes)) > 0)
		genes <- genes[found]
	} else {
		genes <- all_genes[genes]
	}

	# gene with top n norm for each cell
	genes_max <- if (n_top == 1L) {
		all_genes[apply(partials_norm, 1L, function(cell) which.max(cell))]
	} else {
		genes_ord <- t(apply(partials_norm, 1L, function(cell) {
			cell <- setNames(cell, all_genes)
			names(cell[order(cell, decreasing = TRUE)[seq_len(n_top)]])
		}))
		as.vector(genes_ord)
	}
	counts <- as.data.frame(table(genes_max), stringsAsFactors = FALSE)
	counts <- counts[order(counts$Freq, decreasing = TRUE), ]
	counts_valid <- counts[counts$genes_max %in% genes, ]
	gene_ids <- counts_valid$genes_max
	scores <- counts_valid$Freq / sum(counts_valid$Freq)
	names(scores) <- gene_ids
	
	num_top <- min(5L, length(gene_ids))
	top_n_cell_text <- apply(partials_norm, 1L, function(cell) {
		idxs <- head(order(cell, decreasing = TRUE), num_top)
		names <- all_genes[idxs]
		txt <- sprintf('%s. %s (%.3f)', seq_len(num_top), names, cell[idxs])
		paste(txt, collapse = '\n')
	})
	
	# Plot a single map with cells coloured by gene which has 
	# the highest differential norm of all genes considered.
	
	# matrix cells by n_top. might contain NAs later
	max_genes <-
		if (n_top > 1L) genes_ord
		else matrix(gene_ids[apply(partials_norm, 1L, which.max)], ncol = 1)
	# Label smoothing through graph structure
	for (i in seq_len(iter_smooth)) {
		max_genes <- apply(relevance_map@nn_index, 1, function(idx_neighbors) {
			# treating the matrix as vector yields the first column first, i.e. the closest neighbors ...
			max_genes_nn <- unique(as.vector(max_genes[idx_neighbors, ]))
			max_genes_nn_hist <- vapply(max_genes_nn, function(gene) sum(gene == max_genes[idx_neighbors, ], na.rm = TRUE), integer(1L))
			# ... therefore if e.g. all genes only appear once, the first neighbors stay first
			sorted <- sort(max_genes_nn_hist, decreasing = TRUE)
			names(sorted)[seq_len(n_top)]
		})
		# Make inconsistent `apply` result into a cells by n_top matrix
		max_genes <-
			if (n_top > 1) t(max_genes)
			else matrix(max_genes, ncol = 1)
	}
	# use [1] so in case of an empty row we just get an NA.
	gene_labels <- droplevels(structure(
		apply(max_genes, 1, function(cell) na.omit(match(cell, gene_ids))[1]),
		levels = gene_ids,
		class = 'factor'
	))
	# Try creating a palette of the required length
	if (is.function(pal)) pal <- pal(length(unique(gene_labels)))
	if (limit) gene_labels[as.integer(gene_labels) > length(pal)] <- NA

	# Add more than two DC and return data frame so that user
	# can easily rebuild relevance map on other DC combination than 1 and 2.
	rel_map_data <- cbind(as.data.frame(coords), Gene = gene_labels, TopN = top_n_cell_text)
	
	d1 <- colnames(coords)[[1]]
	d2 <- colnames(coords)[[2]]
	rel_map <- ggplot(rel_map_data, aes_string(x = d1, y = d2, colour = 'Gene', text = 'TopN')) +
		geom_point(alpha = .8) +
		geom_rangeframe(colour = par('col')) +
		scale_color_manual(values = pal, na.value = col_na) +
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
