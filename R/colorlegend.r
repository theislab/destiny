#' Color legend
#' 
#' Creates a color legend for a vector used to color a plot. It will use the current \code{\link[grDevices]{palette}()} or the specified \code{pal} as reference.
#' 
#' When passed a factor or integer vector, it will create a discrete legend, whereas a double vector will result in a continuous bar.
#' 
#' @param col          Vector of factor, integer, or double used to determine the ticks.
#' @param pal          If \code{col} is double, pal is used as a continuous palette, else as categorical one
#' @param log          Use logarithmic scale?
#' @param posx         Left and right borders of the color bar relative to plot area (Vector of length 2; 0-1)
#' @param posy         Bottom and top borders of color bar relative to plot area (Vector of length 2; 0-1)
#' @param main         Legend title
#' @param cex_main     Size of legend title font (default: subtitle font size \code{\link{par}('cex.sub')})
#' @param col_main     Color of legend title (default: subtitle color \code{\link{par}('col.sub')})
#' @param col_lab      Color of tick or category labels (default: axis color \code{\link{par}('col.lab')})
#' @param steps        Number of labels in case of a continuous axis. If 0 or FALSE, draw no ticks
#' @param steps_color  Number of gradient samples in case of continuous axis
#' @param digit        Number of digits for continuous axis labels
#' @param left         logical. If TRUE, invert posx
#' @param ...          Additional parameters for the \link[graphics]{text} call used for labels
#' @param cex.main,col.main,col.lab  For compatibility with \code{\link{par}}
#' 
#' @return This function is called for the side effect of adding a colorbar to a plot and returns nothing/NULL.
#' 
#' @examples
#' color_data <- 1:6
#' par(mar = par('mar') + c(0, 0, 0, 3))
#' plot(sample(6), col = color_data)
#' colorlegend(color_data)
#' 
#' @importFrom graphics par rect segments text
#' @importFrom grDevices colorRampPalette palette
#' @export
colorlegend <- function(
	col, pal = palette(), log = FALSE,
	posx = c(.9, .93), posy = c(.05, .9),
	main = NULL, cex_main = par('cex.sub'),
	col_main = par('col.sub'), col_lab = par('col.lab'),
	steps = 5, steps_color = 100,
	digit = 2, left = FALSE,
	...,
	cex.main = NULL,
	col.main = NULL,
	col.lab = NULL) {
	draw_ticks <- as.logical(steps)
	if (!draw_ticks) steps <- 2L
	if (!is.null(cex.main)) cex_main <- cex.main
	if (!is.null(col.main)) col_main <- col.main
	if (!is.null(col.lab))  col_lab  <- col.lab
	
	zval <-
		if      (is.double(col)) seq(min(col, na.rm = TRUE), max(col, na.rm = TRUE), length.out = steps)
		else if (is.factor(col)) factor(levels(col))
		else                     sort(unique(col))
	
	zval_num <-
		if      (is.integer(zval)) seq_along(zval)
		else if (is.numeric(zval)) zval
		else if (is.factor(zval) || is.character(zval)) seq_along(zval)
		else                       as.integer(zval)
	
	zlim <-
		if (is.double(col)) range(zval_num)
		else c(min(zval_num) - .5, max(zval_num) + .5)
	
	par(new = TRUE)
	omar <- nmar <- par('mar')
	nmar[c(2, 4)] <- 0
	par(mar = nmar)
	
	emptyplot()
	
	pars <- par('usr')
	dx <- pars[[2]] - pars[[1]]
	xmin <- pars[[1]] + posx[[1]] * dx
	xmax <- pars[[1]] + posx[[2]] * dx
	dy <- pars[[4]] - pars[[3]]
	ymin <- pars[[3]] + posy[[1]] * dy
	ymax <- pars[[3]] + posy[[2]] * dy
	
	if (log) {
		zlim <- log10(zlim)
		zval <- log10(zval)
	}
	zmin <- zlim[[1]]
	zmax <- zlim[[2]]
	
	if (is.double(col)) {
		batches <- colorRampPalette(pal)(steps_color)
		Y <- seq(ymin, ymax, length.out = length(batches) + 1)
	} else {
		idx_c <- seq(min(zval_num), max(zval_num))
		idx_c[!(idx_c %in% zval_num)] <- NA
		
		batches <- pal[idx_c]
		Y <- seq(ymin, ymax, length.out = length(idx_c) + 1)
	}
	
	rect(xmin, Y[-length(Y)], xmax, Y[-1], col = batches, border = NA)
	rect(xmin, ymin, xmax, ymax, border = col_lab)
	
	dx <- xmax - xmin
	dy <- ymax - ymin
	if (left) {
		Dx <- -dx
		pos <- 2
		xpos <- xmin + Dx * .5
	}
	else {
		Dx <- +dx
		pos <- 4
		xpos <- xmax + Dx * .5
	}
	
	zval_txt <- if (is.double(col)) formatC(zval, digits = digit, format = 'fg') else zval
	
	Ypos <- ymin + (zval_num - zmin)/(zmax - zmin) * dy
	if (draw_ticks) {
		segments(xmax, Ypos, xpos + Dx * .25, Ypos, col = col_lab)
		text(xpos, Ypos, zval_txt, pos = pos, col = col_lab, ...)
	}
	
	if (!is.null(main)) {
		for (i in length(main):1)
			text(x = mean(c(xmin, xmax)),
					 y = ymax + .05 * (length(main) - i + 1),
					 labels = main[i],
					 adj = c(.5, .5),
					 cex = cex_main,
					 col = col_main)
	}
	par(new = FALSE)
	par(mar = omar)
	invisible()
}
