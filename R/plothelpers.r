#' @importFrom graphics par rect segments text
NULL

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
#' @param main.cex     Size of legend title font
#' @param main.col     Color of legend title
#' @param lab.col      Color of tick or category labels
#' @param steps        Number of labels in case of a continuous axis
#' @param color.steps  Number of gradient samples in case of continuous axis
#' @param digit        Number of digits for continuous axis labels
#' @param left         logical. If TRUE, invert posx
#' @param ...          Additional parameters for the \link[graphics]{text} call used for labels
#' 
#' @return This function is called for the side effect of adding a colorbar to a plot and returns nothing/NULL.
#' 
#' @examples
#' color.data <- 1:6
#' par(mar = par('mar') + c(0, 0, 0, 3))
#' plot(sample(6), col = color.data)
#' colorlegend(color.data)
#' 
#' @export
colorlegend <- function(
	col, pal = palette(), log = FALSE,
	posx = c(0.9, 0.93), posy = c(0.05, 0.9),
	main = NULL, main.cex = 1,
	main.col = par('fg'), lab.col = par('fg'),
	steps = 5, color.steps = 100,
	digit = 2, left = FALSE,
...) {
	if (is.double(col)) {
		zval <- seq(min(col), max(col), length.out = steps)
	} else {
		zval <- sort(unique(col))
	}
	
	if (is.integer(zval))
		zval.num <- seq_along(zval)
	else if (is.numeric(zval))
		zval.num <- zval
	else 
		zval.num <- as.integer(zval)
	
	if (is.double(col)) {
		zlim <- range(zval.num)
	} else {
		zlim <- c(min(zval.num) - .5, max(zval.num) + .5)
	}
	
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
		batches <- colorRampPalette(pal)(color.steps)
		Y <- seq(ymin, ymax, length.out = length(batches) + 1)
	} else {
		c.idx <- seq(min(zval.num), max(zval.num))
		c.idx[!(c.idx %in% zval.num)] <- NA
		
		batches <- pal[c.idx]
		Y <- seq(ymin, ymax, length.out = length(c.idx) + 1)
	}
	
	rect(xmin, Y[-length(Y)], xmax, Y[-1], col = batches, border = NA)
	rect(xmin, ymin, xmax, ymax, border = lab.col)
	
	dx <- xmax - xmin
	dy <- ymax - ymin
	if (left) {
		Dx <- -dx
		pos <- 2
		xpos <- xmin + Dx * 0.5
	}
	else {
		Dx <- +dx
		pos <- 4
		xpos <- xmax + Dx * 0.5
	}
	
	zval.txt <- if (is.double(col)) formatC(zval, digits = digit, format = 'fg') else zval
	
	Ypos <- ymin + (zval.num - zmin)/(zmax - zmin) * dy
	segments(xmax, Ypos, xpos + Dx * 0.25, Ypos, col = lab.col)
	text(xpos, Ypos, zval.txt, pos = pos, col = lab.col, ...)
	
	if (!is.null(main)) {
		for (i in length(main):1)
			text(x = mean(c(xmin, xmax)),
					 y = ymax + 0.05 * (length(main) - i + 1),
					 labels = main[i],
					 adj = c(0.5, 0.5),
					 cex = main.cex,
					 col = main.col)
	}
	par(new = FALSE)
	par(mar = omar)
	invisible()
}

#' Sequential color palette using the cube helix system
#' 
#' Creates a perceptually monotonously decreasing (or increasing) lightness color palette with different tones.
#' 
#' @param n        Number of colors to return (default: 6)
#' @param start    Hue to start helix at (\eqn{\textrm{start} \in [0,3]}, default: 0)
#' @param r        Number of rotations of the helix. Can be negative. (default: 0.4)
#' @param hue      Saturation. 0 means greyscale, 1 fully saturated colors (default: 0.8)
#' @param gamma    Emphasize darker (gamma < 1) or lighter (gamma > 1) colors (default: 1)
#' @param light    Lightest lightness (default: 0.85)
#' @param dark     Darkest lightness (default: 0.15)
#' @param reverse  logical. If TRUE, reverse lightness (default: FALSE)
#' 
#' @return A \code{character} vector of hex colors with length \code{n}
#' 
#' @examples
#' palette(cube.helix())
#' image(matrix(1:6), col = 1:6, pch = 19, axes = FALSE)
#' 
#' cr <- colorRamp(cube.helix(12, r = 3))
#' r <- runif(100)
#' plot(1:100, r, col = rgb(cr(r), max = 255), type = 'b', pch = 20)
#' 
#' @export
cube.helix <- function(n = 6, start = 0, r = .4, hue = .8, gamma = 1, light = .85, dark = .15, reverse = FALSE) {
	M <- matrix(c(-0.14861, -0.29227, 1.97294,
	               1.78277, -0.90649, 0), ncol = 2)
	lambda <- seq(light, dark, length.out = n)
	if (reverse) lambda <- rev(lambda)
	l <- rep(lambda ^ gamma, each = 3)
	phi <- 2 * pi * (start/3 + r * lambda)
	t <- rbind(cos(phi), sin(phi))
	out <- l + hue * l * (1 - l)/2 * (M %*% t)
	out <- pmin(pmax(out, 0), 1)
	out <- apply(out, 2, function(x) rgb(x[[1]], x[[2]], x[[3]]))
	out
}

emptyplot <- function(xlim = c(0, 1), ylim = xlim, asp = 1, frame.plot = FALSE, col = NULL, ...) {
	plot(0, type = "n", xlab = "", ylab = "", asp = asp, axes = FALSE, 
	     frame.plot = frame.plot, xlim = xlim, ylim = ylim, xaxs = "i", 
	     yaxs = "i", ...)
	if (!is.null(col)) 
		rect(xlim[1], ylim[1], xlim[2], ylim[2], col = col)
}

continuous.colors <- function(vals, pal = palette(), limits = NULL) {
	if (is.null(limits))
		limits <- range(vals, na.rm = TRUE)
	
	ramp <- colorRamp(pal)
	scaled <- (vals - limits[[1]]) / diff(limits)
	rgb(ramp(scaled), maxColorValue = 255)
}
