#' @importFrom tidyr gather
#' @importFrom scales percent
#' @importFrom ggplot2 ggplot aes_string stat
#' @importFrom ggplot2 scale_fill_gradient scale_alpha_continuous
#' @importFrom ggplot.multistats stat_summaries_hex
#' @rdname Gene-Relevance-plotting
#' @export
plot_rank_bin <- function(gr, gene, ..., n_top = 10L, low = "#3B99B1", high = "#F5191C", bins = 10L) {
	genes_missing <- setdiff(gene, colnames(gr@partials_norm))
	if (length(genes_missing) > 0) {
		genes_close <- lapply(genes_missing, agrep, colnames(gr@partials_norm))
		stop('Missing genes: ', paste(genes_missing, collapse = ', '), '. ',
				 'Closest available: ', paste(unlist(genes_close), collapse = ', '))
	}
	
	top10 <- function(x) sum(x <= 10) / length(x)
	
	partials <- as.data.frame(t(apply(-gr@partials_norm, 1, rank)[gene, , drop = FALSE]))
	coords <- as.data.frame(gr@coords[, 1:2])
	d <- gather(cbind(partials, coords), 'Gene', 'Rank', -tidyselect::one_of(names(coords)))
	
	ggplot(d, aes_string('DC1', 'DC2', z = 'Rank')) +
		stat_summaries_hex(
			aes_string(fill = 'stat(top10)', alpha  = 'stat(size)'),
			funs = list(top10 = top10, size = 'length', 'median'),
			bins = bins
		) +
		scale_fill_gradient(
			name = substitute(Rank <= n_top, list(n_top = n_top)),
			labels = percent,
			low  = low,
			high = high
		) +
		scale_alpha_continuous(name = '#Cells', trans = 'log10') +
		theme_really_minimal()
}
