emptyplot <- function(xlim = c(0, 1), ylim = xlim, asp = 1, frame.plot = FALSE, col = NULL, ...) {
	plot(0, type = 'n', xlab = '', ylab = '', asp = asp, axes = FALSE, 
	     frame.plot = frame.plot, xlim = xlim, ylim = ylim, xaxs = 'i', 
	     yaxs = 'i', ...)
	if (!is.null(col)) 
		rect(xlim[1], ylim[1], xlim[2], ylim[2], col = col)
}

#' @importFrom grDevices palette
#' @importFrom scales colour_ramp rescale
continuous_colors <- function(vals, pal = palette(), limits = NULL, levels = 100) {
	if (is.function(pal))
		pal <- pal(levels)
	
	if (is.null(limits))
		limits <- range(vals, na.rm = TRUE)
	
	ramp <- colour_ramp(pal, alpha = TRUE)
	
	ramp(rescale(vals, from = limits))
}
