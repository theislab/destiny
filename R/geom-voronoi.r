#' Voronoi tesselation
#' 
#' This stat and geom allows you to display voronoi tesselation as polygons.
#' The computations are based on the \code{\link[deldir]{deldir}()} function.
#' 
#' @section Aesthetics:
#' Understands the following aesthetics.
#' Required aesthetics are in bold:
#' 
#' \itemize{
#' \item \strong{x}
#' \item \strong{y}
#' \item alpha
#' \item color
#' \item fill
#' \item linetype
#' \item size
#' }
#' 
#' @inheritParams ggplot2::geom_polygon
#' 
#' @param eps        A value of epsilon used in testing whether a quantity is zero,
#'                   mainly in the context of whether points are collinear.
#'                   If anomalous errors arise, it is possible that these may averted
#'                   by adjusting the value of eps upward or downward.
#' @param normalize  Should coordinates be normalized prior to calculations.
#'                   If \code{x} and \code{y} are in wildly different ranges
#'                   it can lead to tesselation and triangulation that seems off
#'                   when plotted without \code{\link[ggplot2]{coord_fixed}()}.
#'                   Normalization of coordinates solves this.
#'                   The coordinates are transformed back after calculations.
#' @param expand     How much to expand the convex hull around the points.
#'                   Default: 0.1; \eqn{expand \times max(span(data$x), span(data$y))}
#' 
#' @name geom_voronoi
#' @rdname geom_voronoi
#' 
#' @examples
#' library(ggplot2)
#' ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
#'   geom_voronoi(aes(fill = Species)) +
#'   geom_point()
#' 
#' # Difference of normalize = TRUE
#' ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
#'   geom_voronoi(aes(fill = Species), normalize = TRUE) +
#'   geom_point()
NULL


#' @rdname geom_voronoi
#' @importFrom ggplot2 layer GeomPolygon
#' @inheritParams ggplot2::geom_polygon
#' @export
geom_voronoi <- function(
	mapping = NULL, data = NULL, stat = 'voronoi',
	position = 'identity', na.rm = FALSE, eps = 1e-9, normalize = FALSE, expand = .1,
	show.legend = NA, inherit.aes = TRUE, ...
) layer(
	data = data, mapping = mapping, stat = stat, geom = GeomPolygon,
	position = position, show.legend = show.legend, inherit.aes = inherit.aes,
	params = list(eps = eps, normalize = normalize, expand = expand, na.rm = na.rm, ...)
)


#' @rdname geom_voronoi
#' @importFrom scales rescale
#' @importFrom ggplot2 ggproto Stat
#' @importFrom dplyr %>% select
#' @export
StatVoronoi <- ggproto('StatVoronoi', Stat,
	compute_panel = function(data, scales, eps = 1e-9, normalize = FALSE, expand = .1) {
		if (!requireNamespace('sf', quietly = TRUE) || !requireNamespace('lwgeom', quietly = TRUE)) {
			stop('The packages sf and lwgeom are needed for the voronoi stat and geom.')
		}
		
		data$group <- seq_len(nrow(data))
		if (data %>% select('x', 'y') %>% duplicated() %>% any()) {
			warning('stat_voronoi: dropping duplicated points', call. = FALSE)
		}
		if (normalize) {
			x_range <- range(data$x, na.rm = TRUE, finite = TRUE)
			y_range <- range(data$y, na.rm = TRUE, finite = TRUE)
			data$x <- rescale(data$x, from = x_range)
			data$y <- rescale(data$y, from = y_range)
		}
		margin <-
			if (inherits(expand, 'units')) expand
			else expand * max(diff(range(data$x)), diff(range(data$y)))
		tiles <- to_tile(data, margin)
		data$x <- NULL
		data$y <- NULL
		data <- merge(data, tiles, all = TRUE)
		if (normalize) {
			data$x <- rescale(data$x, to = x_range, from = c(0, 1))
			data$y <- rescale(data$y, to = y_range, from = c(0, 1))
		}
		data
	},
	required_aes = c('x', 'y')
)


#' @importFrom dplyr %>% select bind_rows
to_tile <- function(data, margin) {
	points <- data %>% select('x', 'y') %>% as.matrix() %>% sf::st_multipoint()
	hull <- sf::st_convex_hull(points) %>% sf::st_buffer(margin)
	# st_voronoi returns a GEOMETRYCOLLECTION containing only polygons,
	# because a MULTIPOLYGON cannot have shared corner points.
	polys <- sf::st_voronoi(points, hull) %>% sf::st_collection_extract('POLYGON') %>% lwgeom::st_make_valid()
	ord <- points %>% sf::st_sfc() %>% sf::st_cast('POINT') %>% sf::st_intersects(polys) %>% unlist()
	sf::st_intersection(polys[ord], hull) %>%
		lapply(function(poly) poly %>% as.matrix() %>% as.data.frame() %>% setNames(c('x', 'y'))) %>%
		bind_rows(.id = 'group')
}

# ggplot(dm_scial, aes(DC1, DC2, fill = seq_len(ncol(scial)))) +
#   geom_voronoi() + geom_point(shape = 21, colour = 'transparent') +
#   scale_fill_cube_helix(r = 5, discrete = FALSE, name = 'Point') +
#   theme_minimal() + theme(panel.grid = element_blank(), axis.text = element_blank())
