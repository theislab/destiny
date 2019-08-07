#' @rdname Gene-Relevance-plotting
#' @export
setGeneric('plot_differential_map', function(coords, exprs, ..., genes, dims = 1:2, pal = hcl.colors, faceter = facet_wrap(~ Gene)) {
	standardGeneric('plot_differential_map')
})

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot_differential_map', c('matrix', 'matrix'), function(coords, exprs, ..., genes, dims, pal, faceter) {
	plot_differential_map_impl(gene_relevance(coords, exprs, dims = seq_len(max(dims))), genes = genes, dims = dims, pal = pal, faceter = faceter, ...)
})

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot_differential_map', c('DiffusionMap', 'missing'), function(coords, exprs, ..., genes, dims, pal, faceter) {
	plot_differential_map_impl(gene_relevance(coords, dims = seq_len(max(dims))), genes = genes, dims = dims, pal = pal, faceter = faceter, ...)
})

#' @rdname Gene-Relevance-plotting
#' @export
setMethod('plot_differential_map', c('GeneRelevance', 'missing'), function(coords, exprs, ..., genes, dims, pal, faceter) {
	plot_differential_map_impl(coords, genes = genes, dims = dims, pal = pal, faceter = faceter, ...)
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

	do.call(rbind, lapply(genes, function(g) {
		d_var <- .05  # Fraction of overall dimension variability
		partials <- lapply(seq_len(length(dims)), function(d) {
			dc <- coords[, d]
			partials <- relevance_map@partials[, g, d]
			# Scale magnitude of partial derivates
			delta <- diff(rev(range(dc, na.rm = TRUE)))
			partials / max(abs(partials), na.rm = TRUE) * d_var * delta
		})
		
		cbind(
			as.data.frame(coords),
			Expression = exprs[, g],
			PartialsNorm = partials_norms[, g],
			Cell = if (!is.null(rownames(exprs))) rownames(exprs) else seq_len(nrow(exprs)),
			Gene = factor(g, levels = genes),
			Angle     = atan(partials[[2]]   / partials[[1]]  ),
			Magnitude = sqrt(partials[[1]]^2 + partials[[2]]^2))
	}))
}

#' @importFrom ggplot2 ggplot aes aes_string
#' @importFrom ggplot2 geom_point geom_spoke
#' @importFrom ggplot2 scale_colour_gradientn
#' @importFrom ggplot2 ggtitle facet_wrap
#' @importFrom grid arrow unit
plot_differential_map_impl <- function(relevance_map, ..., genes, dims, pal, faceter) {
	chkDots(...)
	if (is.function(pal)) pal <- pal(12)
	dtm <- differential_map(relevance_map, genes, dims)
	coords <- get_coords(relevance_map, dims)
	gene_names <- if (is.character(genes)) genes else colnames(relevance_map@exprs)[genes]
	
	d1 <- colnames(coords)[[1]]
	d2 <- colnames(coords)[[2]]
	gg <- ggplot(dtm, aes_string(d1, d2)) +
		geom_spoke(
			aes_string(
				angle = 'Angle', radius = 'Magnitude',
				alpha = 'PartialsNorm',
				colour = 'Expression'
			),
			arrow = arrow(length = unit(.01, 'npc'))
		) +
		scale_colour_gradientn(colours = pal) + 
		geom_rangeframe(colour = par('col')) +
		theme_really_minimal()
	
	if (length(genes) > 1) gg + faceter else gg + ggtitle(gene_names)
}
