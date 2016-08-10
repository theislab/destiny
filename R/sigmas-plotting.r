#' @include sigmas.r
NULL

langle <- '\u27E8'
rangle <- '\u27E9'
angles <- list(langle = langle, rangle = rangle)

#' Plot \link{Sigmas} object
#' 
#' @param x              Sigmas object to plot
#' @param col            Vector of bar colors or single color for all bars
#' @param col_highlight  Color for highest bar. Overrides col
#' @param col_line       Color for the line and its axis
#' @param type           Plot type of both lines. Can be a vector of length 2 to specify both separately (default: 'b' aka \dQuote{both lines and points})
#' @param pch            Point identifier for both lines. Can be a vector of length 2 to specify both separately (default: \code{par(pch)} and 4 (a \sQuote{\eqn{\times}}))
#' @param only_dim       logical. If TRUE, only plot the derivative line
#' @param ...            Options passed to the call to plot
#' @param xlab           X label. NULL to use default
#' @param ylab           Either one y label or y labels for both plots. NULL to use both defauts, a NULL in a list of length 2 to use one default.
#' @param main           Title of the plot
#' 
#' @return This method plots a Sigma object to the current device and returns nothing/NULL
#' 
#' @examples
#' data(guo)
#' sigs <- find_sigmas(guo)
#' plot(sigs)
#' 
#' @aliases plot.Sigmas plot,Sigmas,missing-method
#' 
#' @importFrom graphics plot plot.window par text dotchart layout
#' @name plot.Sigmas
#' @export
setMethod('plot', c(x = 'Sigmas', y = 'missing'), function(
	x,
	col = par('fg'),
	col_highlight = '#E41A1C',  # brewer Set1[[1]]
	col_line      = '#999999',  # brewer Set1[[9]]
	type = c('b', 'b'),
	pch = c(par('pch'), 4L),
	only_dim = FALSE,
	...,
	xlab = NULL,
	ylab = NULL,
	main = ''
) {
	sigmas <- x
	if (is.null(optimal_sigma(sigmas))) {
		plot.new()
		plot.window(c(-1,1), c(-1,1))
		text(0, 0, expression('local ' * sigma), .5)
		return()
	}
	if (is.null(sigmas@log_sigmas)) {
		return(dotchart(optimal_sigma(sigmas), labels = 'optimal_sigma', color = col_highlight, pch = pch, main = 'optimal sigma given directly'))
	}
	
	steps <- length(sigmas@log_sigmas)
	
	colors <- if (length(col) == 1L) rep(col, steps - 1L) else col
	colors[[sigmas@optimal_idx]] <- col_highlight
	
	#prepare parameters and reset them at the end
	
	mar <- par('mar')
	
	if (!only_dim) {
		layout(matrix(1:2, 2L))
		par_old <- par(mar = c(0, mar[2:4]), ..., no.readonly = TRUE)
		
		on.exit({
			layout(matrix(1))
			par(par_old)
		})
	}
	
	if (is.null(xlab)) xlab <- expression(log[10](sigma))
	
	if (is.null(ylab)) ylab <- list(NULL, NULL)
	if (length(ylab) == 1L) ylab <- list(ylab, NULL)
	if (is.null(ylab[[1L]]))
		ylab[[1L]] <- substitute(paste(langle, d, rangle), angles)
	if (is.null(ylab[[2L]]))
		ylab[[2L]] <- substitute(paste(langle, log[10](Z(x)), rangle), angles)
	
	#first plot
	
	x <- sigmas@log_sigmas[-1L] - diff(sigmas@log_sigmas) / 2
	xlim <- range(sigmas@log_sigmas)
	
	ysteps <- pretty(sigmas@dim_norms)
	ymin <- min(ysteps)
	ymax <- max(ysteps) + diff(ysteps[1:2])/4
	
	plot(x, sigmas@dim_norms,
			 xlim = xlim, ylim = c(ymin, ymax),
			 xaxt = if (only_dim) 's' else 'n', xlab = if (only_dim) xlab else '',
			 ylab = ylab[[1L]],
			 col = colors, type = type[[1L]], pch = pch[[1L]],
			 main = main,
			 ...)
	
	text(log10(optimal_sigma(sigmas)), sigmas@dim_norms[[sigmas@optimal_idx]],
			 substitute(sigma == rsig, list(rsig = round(optimal_sigma(sigmas), 1))), pos = 3)
	
	#second (overlay) plot
	
	if (!only_dim) {
		par(mar = c(mar[1:2], 0, mar[[4L]]), ...)
		
		col_lines <- rep(col_line, length(sigmas@log_sigmas))
		col_lines[c(sigmas@optimal_idx, sigmas@optimal_idx + 1L)] <- col_highlight
		
		t2 <- if (length(type) > 1L) type[[2L]] else type
		p2 <- if (length(pch)  > 1L)  pch[[2L]] else pch
		plot(sigmas@log_sigmas, sigmas@avrd_norms,
				 xlab = xlab, xlim = xlim,
				 ylab = ylab[[2L]],
				 type = t2, pch = p2, col = col_lines,
				 ...)
	}
	invisible()
})
