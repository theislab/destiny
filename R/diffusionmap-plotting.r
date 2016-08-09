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
#' @param new.dcs      An optional matrix also containing the rows specified with \code{y} and plotted. (default: no more points)
#' @param col          Single color string or vector of discrete or categoric values to be mapped to colors.
#'                     E.g. a column of the data matrix used for creation of the diffusion map. (default: \link{par}\code{('fg')})
#' @param col.by       Specify a \code{dataset(x)} or \code{phenoData(dataset(x))} column to use as color
#' @param col.limits   If \code{col} is a continuous (=double) vector, this can be overridden to map the color range differently than from min to max (e.g. specify \code{c(0, 1)})
#' @param col.new      If \code{new.dcs} is given, it will take on this color. (default: red)
#' @param pal          Palette used to map the \code{col} vector to colors (default: \link{palette}\code{()})
#' @param ...          Parameters passed to \link{plot}, \link[scatterplot3d]{scatterplot3d}, or \link[rgl]{plot3d} (if \code{interactive == TRUE})
#' @param mar          Bottom, left, top, and right margins (default: \code{par(mar)})
#' @param ticks        logical. If TRUE, show axis ticks (default: FALSE)
#' @param axes         logical. If TRUE, draw plot axes (default: Only if \code{tick.marks} is TRUE)
#' @param box          logical. If TRUE, draw plot frame (default: TRUE or the same as \code{axes} if specified)
#' @param legend.main  Title of legend. (default: nothing unless col.by is given)
#' @param legend.opts  Other \link{colorlegend} options (default: empty list)
#' @param interactive  Use \link[rgl]{plot3d} to plot instead of \link[scatterplot3d]{scatterplot3d}?
#' @param draw.legend  logical. If TRUE, draw color legend (default: TRUE if \code{col} is given and a vector to be mapped)
#' @param consec.col   If \code{col} or \code{col.by} refers to an integer column, with gaps (e.g. \code{c(5,0,0,3)}) use the palette color consecutively (e.g. \code{c(3,1,1,2)})
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
#' @importFrom grDevices palette
#' @importFrom scatterplot3d scatterplot3d
#' 
#' @name plot.DiffusionMap
#' @export
plot.DiffusionMap <- function(
	x, dims,
	new.dcs = NULL,
	col = NULL, col.by = NULL, col.limits = NULL,
	col.new = 'red',
	pal = palette(),
	...,
	mar = NULL,
	ticks = FALSE,
	axes = TRUE,
	box = FALSE,
	legend.main = col.by, legend.opts = list(),
	interactive = FALSE,
	draw.legend = !is.null(col) && length(col) > 1 && !is.character(col),
	consec.col = TRUE
) {
	dif <- x
	
	if (interactive) {
		if (!requireNamespace('rgl', quietly = TRUE))
			stop(sprintf('The package %s is required for interactive plots', sQuote('rgl')))
		if (length(dims) != 3L)
			stop('Only 3d plots can be made interactive')
	}
	
	flip <- dims < 0
	dims[flip] <- -dims[flip]
	
	if (!is.null(col) && !is.null(col.by)) stop('Only specify one of col or col.by')
	
	default.col <- is.null(col) && is.null(col.by)
	if (default.col) col <- par('col')
	
	#get col from data
	if (!is.null(col.by))
		col <- extract.col(dataset(dif), col.by)
	
	# extend margin for legend
	old.mar <- par('mar')
	if (is.null(mar)) {
		mar <- old.mar
		if (draw.legend) {
			mar[[4]] <- mar[[4]] + 5
		}
	}
	par(mar = mar)
	
	# reduce axis label position if it is the default and no tick marks are given
	if (all(par('mgp') == c(3, 1, 0)) && !ticks) {
		on.exit(par(mgp = c(3, 1, 0)))
		par(mgp = c(1, 0, 0))
	}
	
	#make consecutive the colors for the color legend
	if (is.integer(col) && consec.col) {
		# c(5,0,0,3) -> c(3,1,1,2)
		col <- factor(col)
	}
	
	#colors to assign to plot rows
	col.plot <- col
	if (is.double(col))
		col.plot <- continuous.colors(col, pal, col.limits)
	else if (is.factor(col.plot))
		col.plot <- as.integer(col.plot)
	
	#limit pal to number of existing colors to maybe attach col.new
	pal.length <- if (is.factor(col)) length(levels(col)) else length(unique(col.plot))
	if (is.function(pal)) {
		# pal is a colorRampPalette-type function
		pal <- pal(pal.length)
	} else {
		# pal is a vector
		pal.length <- min(length(pal), pal.length)
		pal <- pal[seq_len(pal.length)]
	}
	
	#set col.plot to color strings
	if (is.integer(col.plot)) {
		wrapped.idx <- ((col.plot - 1L) %% pal.length) + 1L
		col.plot <- pal[wrapped.idx]
	}
	
	#attach the new.dcs and col.new parameters to data and colors
	point.data <- eigenvectors(dif)[, dims]
	point.data[, flip] <- -point.data[, flip]
	if (!is.null(new.dcs)) {
		point.data <- rbind(point.data, as.matrix(new.dcs[, dims]))
		if (!is.character(col.new)) {
			wrapped.idx <- ((col.plot - 1L) %% pal.length) + 1L
			col.new <- pal[wrapped.idx]
		}
		col.plot <- c(rep_len(col.plot, nrow(point.data) - nrow(new.dcs)),
									rep_len(col.new, nrow(new.dcs)))
		#colorlegend
		if (!is.double(col)) {
			col <- factor(c(as.character(col), 'proj'))
		}
		pal <- c(pal, col.new)
	}
	
	if (length(dims) == 2) {
		p <- NULL
		
		plot(point.data, ..., col = col.plot, axes = FALSE, frame.plot = box)
		if (ticks) {
			r1 <- NULL
			r2 <- NULL
			tl <- 1
		} else {
			r1 <- range(point.data[, 1L])
			r2 <- range(point.data[, 2L])
			tl <- 0
		}
		al <- if (axes && !box) 1 else 0
		axis(1, r1, labels = ticks, lwd = al, lwd.ticks = tl)
		axis(2, r2, labels = ticks, lwd = al, lwd.ticks = tl)
	} else if (length(dims) == 3L) {
		if (interactive) {
			p <- rgl::plot3d(point.data, ..., col = col.plot, axes = FALSE, box = FALSE)
			if (axes || ticks) {
				axtype = if (axes) 'lines' else 'cull'
				nticks = if (ticks) 5 else 0
				rgl::bbox3d(xlen = nticks, ylen = nticks, zlen = nticks, front = axtype, back = axtype)
			}
			if (box) rgl::box3d()
		} else {
			p <- scatterplot3d(
				point.data, ..., color = col.plot, mar = mar,
				axis = axes || box || ticks, lty.axis = if (axes || box) 'solid' else 'blank',
				box = box, tick.marks = ticks)
		}
	} else stop(sprintf('dims is of wrong length (%s): Can only handle 2 or 3 dimensions', dims))
	
	if (draw.legend) {
		legend.col <- if (is.double(col) && !is.null(col.limits)) col.limits else col
		args <- c(list(legend.col, pal = pal, main = legend.main), legend.opts)
		if (interactive) {
			rgl::bgplot3d({
				plot.new()
				do.call(colorlegend, args)
			})
		} else {
			do.call(colorlegend, args)
		}
	}
	
	par(mar = old.mar)
	invisible(p)
}

#' @importFrom Biobase varLabels exprs
extract.col <- function(annot.data, col.by) tryCatch({
	if (inherits(annot.data, 'ExpressionSet')) {
		if (col.by %in% varLabels(annot.data))
			annot.data[[col.by]]
		else
			exprs(annot.data)[col.by, ]
	} else {
		annot.data[, col.by]
	}
}, error = function(e) stop(sprintf('Invalid `col.by`: No column, annotation, or feature found with name %s', dQuote(col.by))))

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
