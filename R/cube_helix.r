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
#' @param discrete If TRUE, return a discrete scale, if FALSE a continuous one (default: TRUE)
#' @param guide    Type of scale guide to use. See \code{\link[ggplot2]{guides}}
#' @param ...      parameters passed to \code{\link[ggplot2]{discrete_scale}} or \code{\link[ggplot2]{continuous_scale}}
#' 
#' @return A \code{character} vector of hex colors with length \code{n}
#' 
#' @examples
#' palette(cube_helix())
#' image(matrix(1:6), col = 1:6, pch = 19, axes = FALSE)
#' 
#' cr <- scales::colour_ramp(cube_helix(12, r = 3))
#' r <- runif(100)
#' plot(1:100, r, col = cr(r), type = 'b', pch = 20)
#' 
#' @importFrom grDevices rgb
#' @export
cube_helix <- function(n = 6, start = 0, r = .4, hue = .8, gamma = 1, light = .85, dark = .15, reverse = FALSE) {
	M <- matrix(c(-.14861, -.29227, 1.97294,
								1.78277, -.90649, 0), ncol = 2)
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

scale_cube_helix <- function(aesthetics, start, r, hue, gamma, light, dark, reverse, discrete, guide, ...) {
	if (!requireNamespace('ggplot2', quietly = TRUE))
		stop('scale_', aesthetics, '_cube_helix needs (and is only useful for) the ggplot2 package')
	
	f <- function(n) cube_helix(n, start, r, hue, gamma, light, dark, reverse)
	
	if (discrete) {
		ggplot2::discrete_scale('colour', 'cube_helix', f, ..., guide = guide)
	} else {
		ggplot2::continuous_scale('colour', 'cube_helix', scales::gradient_n_pal(f(100)), guide = guide, ...)
	}
}

#' @name cube_helix
#' @export
scale_colour_cube_helix <- function(...,
	start = 0, r = .4, hue = .8, gamma = 1, light = .85, dark = .15, reverse = FALSE,
	discrete = TRUE, guide = if (discrete) 'legend' else 'colourbar'
) {
	scale_cube_helix('colour', start, r, hue, gamma, light, dark, reverse, discrete, guide, ...)
}
#' @name cube_helix
#' @export
scale_color_cube_helix <- scale_colour_cube_helix

#' @name cube_helix
#' @export
scale_fill_cube_helix <- function(...,
	start = 0, r = .4, hue = .8, gamma = 1, light = .85, dark = .15, reverse = FALSE,
	discrete = TRUE, guide = if (discrete) 'legend' else 'colourbar'
) {
	scale_cube_helix('fill', start, r, hue, gamma, light, dark, reverse, discrete, guide, ...)
}
