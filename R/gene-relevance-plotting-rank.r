#' @rdname Gene-Relevance-plotting
#' @export
setGeneric('plot_gene_relevance_rank', function(coords, exprs, ..., genes, dims = 1:2, n_top = 10L, pal = c("#3B99B1", "#F5191C"), bins = 10L, faceter = facet_wrap(~ Gene)) {
	standardGeneric('plot_gene_relevance_rank')
})

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot_gene_relevance_rank', c('matrix', 'matrix'), function(coords, exprs, ..., genes, dims, n_top, pal, bins, faceter) {
	plot_gene_relevance_rank_impl(gene_relevance(coords, exprs, dims = seq_len(max(dims))), genes = genes, dims = dims, n_top = n_top, pal = pal, bins = bins, faceter = faceter, ...)
})

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot_gene_relevance_rank', c('DiffusionMap', 'missing'), function(coords, exprs, ..., genes, dims, n_top, pal, bins, faceter) {
	plot_gene_relevance_rank_impl(gene_relevance(coords, dims = seq_len(max(dims))), genes = genes, dims = dims, n_top = n_top, pal = pal, bins = bins, faceter = faceter, ...)
})

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot_gene_relevance_rank', c('GeneRelevance', 'missing'), function(coords, exprs, ..., genes, dims, n_top, pal, bins, faceter) {
	plot_gene_relevance_rank_impl(coords, genes = genes, dims = dims, n_top = n_top, pal = pal, bins = bins, faceter = faceter, ...)
})

#' @importFrom tidyr gather
#' @importFrom scales percent
#' @importFrom ggplot2 ggplot aes_string stat
#' @importFrom ggplot2 scale_fill_gradientn scale_alpha_continuous
#' @importFrom ggplot.multistats stat_summaries_hex
plot_gene_relevance_rank_impl <- function(relevance_map, ..., genes, dims, n_top, pal, bins, faceter) {
	stopifparams(...)
	if (is.function(pal)) pal <- pal(12)
	coords <- get_coords(relevance_map, dims)
	gene_names <- if (is.character(genes)) genes else colnames(relevance_map@exprs)[genes]
	
	genes_missing <- setdiff(genes, colnames(relevance_map@partials_norm))
	if (length(genes_missing) > 0) {
		genes_close <- lapply(genes_missing, agrep, colnames(relevance_map@partials_norm), value = TRUE)
		stop('Missing genes: ', paste(genes_missing, collapse = ', '), '. ',
				 'Closest available: ', paste(unlist(genes_close), collapse = ', '))
	}
	
	top10 <- function(x) sum(x <= 10) / length(x)
	
	partials <- as.data.frame(t(apply(-relevance_map@partials_norm, 1, rank)[genes, , drop = FALSE]))
	d <- gather(cbind(partials, as.data.frame(coords)), 'Gene', 'Rank', !!gene_names)
	
	d1 <- colnames(coords)[[1]]
	d2 <- colnames(coords)[[2]]
	gg <- ggplot(d, aes_string(d1, d2, z = 'Rank')) +
		stat_summaries_hex(
			aes_string(fill = 'stat(top10)', alpha  = 'stat(size)'),
			funs = list(top10 = top10, size = 'length', 'median'),
			bins = bins
		) +
		scale_fill_gradientn(
			name = substitute(Rank <= n_top, list(n_top = n_top)),
			labels = percent,
			colours = pal
		) +
		scale_alpha_continuous(name = '#Cells', trans = 'log10') +
		theme_really_minimal()
	
	if (length(genes) > 1) gg + faceter else gg + ggtitle(gene_names)
}
