#' @include diffusionmap.r
#' @include plothelpers.r
NULL

#' 3D or 2D plot of diffusion map
#' 
#' If you want to plot the eigenvalues, simply \code{plot(eigenvalues(dm)[start:end], ...)}
#' 
#' If you specify negative numbers as diffusion components (e.g. \code{plot(dm, c(-1,2))}), then the corresponding components will be flipped.
#' 
#' @param x            A \link{DiffusionMap}
#' @param dims,y       Diffusion components (eigenvectors) to plot (default: first three components; 1:3)
#' @param new_dcs      An optional matrix also containing the rows specified with \code{y} and plotted. (default: no more points)
#' @param new_data     A data set in the same format as \code{x} that is used to create \code{new_dcs <- \link{dm_predict}(dif, new_data)}
#' @param col          Single color string or vector of discrete or categoric values to be mapped to colors.
#'                     E.g. a column of the data matrix used for creation of the diffusion map. (default: \code{\link[igraph]{cluster_louvain}})
#' @param col_by       Specify a \code{dataset(x)} or \code{phenoData(dataset(x))} column to use as color
#' @param col_limits   If \code{col} is a continuous (=double) vector, this can be overridden to map the color range differently than from min to max (e.g. specify \code{c(0, 1)})
#' @param col_new      If \code{new_dcs} is given, it will take on this color. A vector is also possible. (default: red)
#' @param pal          Palette used to map the \code{col} vector to colors. (default: use \code{\link{cube_helix}} for continuous and \code{\link{palette}()} for discrete data)
#' @param pal_new      Palette used to map the \code{col_new} vector to colors. (default: see \code{pal} argument)
#' @param ...          Parameters passed to \link{plot}, \link[scatterplot3d]{scatterplot3d}, or \link[rgl]{plot3d} (if \code{interactive == TRUE})
#' @param ticks        logical. If TRUE, show axis ticks (default: FALSE)
#' @param axes         logical. If TRUE, draw plot axes (default: Only if \code{ticks} is TRUE)
#' @param box          logical. If TRUE, draw plot frame (default: TRUE or the same as \code{axes} if specified)
#' @param legend_main  Title of legend. (default: nothing unless \code{col_by} is given)
#' @param legend_opts  Other \link{colorlegend} options (default: empty list)
#' @param interactive  Use \link[rgl]{plot3d} to plot instead of \link[scatterplot3d]{scatterplot3d}?
#' @param draw_legend  logical. If TRUE, draw color legend (default: TRUE if \code{col_by} is given or \code{col} is given and a vector to be mapped)
#' @param consec_col   If \code{col} or \code{col_by} refers to an integer column, with gaps (e.g. \code{c(5,0,0,3)}) use the palette color consecutively (e.g. \code{c(3,1,1,2)})
#' @param col_na       Color for \code{NA} in the data. specify \code{NA} to hide.
#' @param plot_more    Function that will be called while the plot margins are temporarily changed
#'                     (its \code{p} argument is the rgl or scatterplot3d instance or NULL,
#'                     its \code{rescale} argument is \code{NULL}, a \code{list(from = c(a, b), to = c(c, d))}),
#'                     or an array of shape \eqn{from|to \times dims \times min|max}, i.e. \eqn{2 \times length(dims) \times 2}.
#'                     In case of 2d plotting, it should take and return a ggplot2 object.
#' 
#' @return The return value of the underlying call is returned, i.e. a scatterplot3d or rgl object.
#' 
#' @examples
#' data(guo)
#' plot(DiffusionMap(guo))
#' 
#' @aliases plot.DiffusionMap plot,DiffusionMap,missing-method plot,DiffusionMap,numeric-method
#' 
#' @importFrom graphics par axis plot plot.new
#' @importFrom stats setNames
#' @importFrom grDevices palette
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom ggplot2 ggplot aes aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme theme_minimal element_blank element_line element_text element_rect
#' @importFrom ggplot2 scale_fill_identity scale_fill_manual scale_fill_gradientn scale_fill_identity
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 guide_colourbar guide_legend
#' @importFrom ggthemes geom_rangeframe extended_range_breaks
#' 
#' @name plot.DiffusionMap
#' @export
plot.DiffusionMap <- function(
	x, dims = 1:3,
	new_dcs = if (!is.null(new_data)) dm_predict(x, new_data),
	new_data = NULL,
	col = NULL, col_by = NULL, col_limits = NULL,
	col_new = 'red',
	pal = NULL, pal_new = NULL,
	...,
	ticks = FALSE,
	axes = TRUE,
	box = FALSE,
	legend_main = col_by, legend_opts = list(),
	interactive = FALSE,
	draw_legend = !is.null(col_by) || (length(col) > 1 && !is.character(col)),
	consec_col = TRUE, col_na = 'grey',
	plot_more = function(p, ..., rescale = NULL) p
) {
	dif <- x
	is_projection <- !is.null(new_dcs) && is.character(col_new) && length(col_new) == 1L
	
	if (interactive) {
		if (!requireNamespace('rgl', quietly = TRUE))
			stop(sprintf('The package %s is required for interactive plots', sQuote('rgl')))
		if (length(dims) != 3L)
			stop('Only 3d plots can be made interactive')
	}
	
	if (!is.null(col) && !is.null(col_by)) stop('Only specify one of col or col_by')
	if (!is.null(col_by)) {
		col <- dataset_get_feature(dataset(dif), col_by)
	} else if (is.null(col)) {
		col <- get_louvain_clusters(dif@transitions)
		col_by <- legend_main <- 'Louvain'
	}
	continuous <- is.double(col)
	if (is_projection) {
		projection_guide <- setNames(c(col, col_new), c(paste(legend_main, col), rep_len('new', length(col_new))))
		legend_main <- 'Projection'
	}
	col_legend <- if (continuous && !is.null(col_limits)) col_limits else col
	
	# use a fitting default palette
	if (is.null(pal)) {
		pal <- if (is.double(col)) cube_helix
		else palette()
	}
	
	# make consecutive the colors for the color legend
	if (is.integer(col) && consec_col) {
		# c(5,0,0,3) -> c(3,1,1,2)
		col <- factor(col)
	}
	
	point_data <- cbind(
		as.data.frame(flipped_dcs(eigenvectors(dif), dims)),
		Colour     = col,
		ColourExpl = get_explicit_col(col, pal, col_na, col_limits),
		Projection = factor(rep('old', nrow(eigenvectors(dif))), c('old', 'new')))
	rm(col)
	
	if (!is.null(new_dcs)) {
		point_data <- rbind(point_data, cbind(
			as.data.frame(flipped_dcs(new_dcs, dims)),
			Colour     = col_new,  #TODO
			ColourExpl = get_explicit_col(col_new, pal_new, col_na, col_limits),
			Projection = 'new'
		))
		col_legend
	}
	
	col_lvls <- na.omit(as.character(unique(point_data$Colour)))
	is_one_colour <- length(col_lvls) == 1L
	
	if (length(dims) == 2) {
		d1 <- names(point_data)[[1L]]
		d2 <- names(point_data)[[2L]]
		
		p <- ggplot(point_data, aes_string(d1, d2)) +
			theme_minimal() + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
		
		# TODO: this logic might be off
		if (is_projection || !is_one_colour) p <- p +
			geom_point(
				aes_string(fill = if (continuous || is_projection || !is.null(col_by)) 'Colour' else 'ColourExpl'),
				colour = I('transparent'),
				shape  = I(21))
		
		if (!is_one_colour) p <- p +
			if (is_projection)   scale_fill_identity (name = legend_main, guide = 'legend', labels = names(projection_guide), breaks = projection_guide, na.value = col_na)
			else if (continuous) scale_fill_gradientn(name = legend_main, colours = if (is.function(pal)) pal(100) else pal, na.value = col_na)
			else                 scale_fill_manual   (name = legend_main, values  = c(if (is.function(pal)) pal(length(col_lvls)) else pal[seq_along(col_lvls)], NA), labels = c(col_lvls, NA), na.value = col_na)
		if (box)   p <- p + theme(panel.border = element_rect(fill = NA), axis.title.x = element_text(), axis.title.y = element_text())
		if (ticks) p <- p + theme(axis.ticks = element_line(), axis.text.x  = element_text(), axis.text.y  = element_text())
		if (axes)  p <- p + geom_rangeframe(colour = par('col'))
		if (ticks && axes && !box) p <- p + 
			scale_x_continuous(breaks = extended_range_breaks()(point_data[[1L]])) +
			scale_y_continuous(breaks = extended_range_breaks()(point_data[[2L]]))
		p <- plot_more(p, rescale = NULL)
	} else if (length(dims) == 3L) {
		if (interactive) {
			p <- rgl::plot3d(point_data, ..., col = point_data$ColourExpl, axes = FALSE, box = FALSE)
			if (axes || ticks) {
				axtype = if (axes) 'lines' else 'cull'
				nticks = if (ticks) 5 else 0
				rgl::bbox3d(xlen = nticks, ylen = nticks, zlen = nticks, front = axtype, back = axtype)
			}
			if (box) rgl::box3d()
			plot_more(p, rescale = NULL)
		} else {
			rescale <- NULL
			if (!ticks) {
				rescale <- array(NA, c(2L, length(dims), 2L), list(c('from', 'to'), as.character(dims), c('min', 'max')))
				for (d in seq_along(dims)) {  # -> scatterplot3d's pretty() should not mess things up
					r <- range(point_data[, d])
					point_data[, d] <- scales::rescale(point_data[, d], c(0, 1), r)
					rescale['from', d, ] <- r
					rescale['to', d, ] <- c(0, 1)
				}
			}
			
			mar <- list(...)$mar
			if (is.null(mar)) mar <- par('mar')
			old_mar <- mar; on.exit(par(mar = old_mar))
			if (draw_legend) mar[[4]] <- mar[[4]] + 5
			p <- scatterplot3d(
				point_data[, 1:3], ..., color = point_data$ColourExpl, mar = mar,
				axis = axes || box || ticks, lty.axis = if (axes || box) 'solid' else 'blank',
				box = box, tick.marks = ticks)
			rm(mar)
			plot_more(p, rescale = rescale)
			
			if (draw_legend) {
				args <- c(list(col_legend, pal = pal, main = legend_main), legend_opts)
				if (interactive) {
					rgl::bgplot3d({
						plot.new()
						do.call(colorlegend, args)
					})
				} else {
					do.call(colorlegend, args)
				}
			}
		}
	} else stop(sprintf('dims is of wrong length (%s): Can only handle 2 or 3 dimensions', dims))
	
	if (length(dims) == 2) p else invisible(p)
}


get_explicit_col <- function(col, pal, col_na, col_limits) {
	# if nothing is given, return one colour
	if (is.null(col)) return(par('col'))
	
	# if we have continuous colour, we are done.
	if (is.double(col))
		return(continuous_colors(col, pal, col_limits))
	
	# get palette length and convert col to consecutive integers
	length_pal <-
		if (is.factor(col))
			length(levels(col))
		else if (is.integer(col))
			length(unique(col))
		else stopifnot(is.character(col))
	if (is.factor(col))
		col <- as.integer(col)
	
	# map integers to strings if necessary
	if (is.integer(col)) {
		if (is.function(pal)) {
			# pal is a colorRampPalette-type function
			pal <- pal(length_pal)
		} else {
			# pal is a vector
			length_pal <- min(length(pal), length_pal)
			pal <- pal[seq_len(length_pal)]
		}
		
		idx_wrapped <- ((col - 1L) %% length_pal) + 1L
		col <- pal[idx_wrapped]
		col[is.na(col)] <- col_na
	}
	
	# if the color wasn’t numeric, use as is
	col
}

# test:
# layout(matrix(1:8, 2))
# mapply(function(t, a, b) plot(dif, ticks = t, axes = a, box = b, main = sprintf('t=%s a=%s b=%s', t, a, b)),
#        c(T,T,T,T,F,F,F,F), c(T,F,T,F,T,F,T,F), c(T,T,F,F,T,T,F,F))

#' @name plot.DiffusionMap
#' @export
setMethod('plot', c(x = 'DiffusionMap', y = 'numeric'), function(x, y, ...) plot.DiffusionMap(x, y, ...))

#' @name plot.DiffusionMap
#' @export
setMethod('plot', c(x = 'DiffusionMap', y = 'missing'), function(x, y, ...) plot(x, seq_len(min(3L, ncol(eigenvectors(x)))), ...))
