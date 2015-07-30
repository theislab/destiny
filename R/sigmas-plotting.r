#' @include sigmas.r
NULL

langle <- '\u27E8'
rangle <- '\u27E9'
angles <- list(langle = langle, rangle = rangle)

#' Plot \link{Sigmas} object
#' 
#' @param x              Sigmas object to plot
#' @param col            Vector of bar colors or single color for all bars
#' @param highlight.col  Color for highest bar. Overrides col
#' @param line.col       Color for the line and its axis
#' @param type           Plot type of both lines. Can be a vector of length 2 to specify both separately (default: 'b' aka \dQuote{both lines and points})
#' @param pch            Point identifier for both lines. Can be a vector of length 2 to specify both separately (default: \code{par(pch)} and 4 (a \sQuote{\eqn{\times}}))
#' @param only.dim       logical. If TRUE, only plot the derivative line
#' @param ...            Options passed to the call to plot
#' @param xlab           X label. NULL to use default
#' @param ylab           Either one y label or y labels for both plots. NULL to use both defauts, a NULL in a list of length 2 to use one default.
#' @param main           Title of the plot
#' 
#' @return This method plots a Sigma object to the current device and returns nothing/NULL
#' 
#' @examples
#' data(guo)
#' sigs <- find.sigmas(guo)
#' plot(sigs)
#' 
#' @aliases plot.Sigmas plot,Sigmas,missing-method
#' @name plot.Sigmas
#' @importFrom graphics plot
#' @export
setMethod('plot', c(x = 'Sigmas', y = 'missing'), function(
	x,
	col = par('fg'),
	highlight.col = '#E41A1C',  # brewer Set1[[1]]
	line.col      = '#999999',  # brewer Set1[[9]]
	type = c('b', 'b'),
	pch = c(par('pch'), 4L),
	only.dim = FALSE,
	...,
	xlab = NULL,
	ylab = NULL,
	main = ''
) {
	sigmas <- x
	if (is.null(sigmas@log.sigmas)) {
		return(dotchart(sigmas@optimal.sigma, labels = 'optimal.sigma', color = highlight.col, pch = pch, main = 'optimal sigma given directly'))
	}
	
	steps <- length(sigmas@log.sigmas)
	
	colors <- if (length(col) == 1L) rep(col, steps - 1L) else col
	colors[[sigmas@optimal.idx]] <- highlight.col
	
	#prepare parameters and reset them at the end
	
	mar <- par('mar')
	
	if (!only.dim) {
		layout(matrix(1:2, 2L))
		old.par <- par(mar = c(0, mar[2:4]), ..., no.readonly = TRUE)
		
		on.exit({
			layout(matrix(1))
			par(old.par)
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
	
	x <- sigmas@log.sigmas[-1L] - diff(sigmas@log.sigmas) / 2
	xlim <- range(sigmas@log.sigmas)
	
	ysteps <- pretty(sigmas@dim.norms)
	ymin <- min(ysteps)
	ymax <- max(ysteps) + diff(ysteps[1:2])/4
	
	plot(x, sigmas@dim.norms,
			 xlim = xlim, ylim = c(ymin, ymax),
			 xaxt = if (only.dim) 's' else 'n', xlab = if (only.dim) xlab else '',
			 ylab = ylab[[1L]],
			 col = colors, type = type[[1L]], pch = pch[[1L]],
			 main = main,
			 ...)
	
	text(log10(sigmas@optimal.sigma), sigmas@dim.norms[[sigmas@optimal.idx]],
			 substitute(sigma == rsig, list(rsig = round(sigmas@optimal.sigma, 1))), pos = 3)
	
	#second (overlay) plot
	
	if (!only.dim) {
		par(mar = c(mar[1:2], 0, mar[[4L]]), ...)
		
		l.cols <- rep(line.col, length(sigmas@log.sigmas))
		l.cols[c(sigmas@optimal.idx, sigmas@optimal.idx + 1L)] <- highlight.col
		
		t2 <- if (length(type) > 1L) type[[2L]] else type
		p2 <- if (length(pch)  > 1L)  pch[[2L]] else pch
		plot(sigmas@log.sigmas, sigmas@avrd.norms,
				 xlab = xlab, xlim = xlim,
				 ylab = ylab[[2L]],
				 type = t2, pch = p2, col = l.cols,
				 ...)
	}
	invisible()
})
