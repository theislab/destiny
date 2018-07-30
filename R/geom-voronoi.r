#' Voronoi tesselation
#'
#' This stat and geom allows you to display voronoi tesselation as polygons.
#' The computations are based on the \code{\link[deldir]{deldir}()} function.
#'
#' @section Aesthetics:
#' Understands the following aesthetics.
#' Required aesthetics are in bold:
#'
#' - **x**
#' - **y**
#' - alpha
#' - color
#' - fill
#' - linetype
#' - size
#'
#' @inheritParams ggplot2::geom_polygon
#'
#' @param bound The bounding rectangle for the tesselation. Defaults to
#' \code{NULL} which creates a rectangle expanded 10% in all directions. If
#' supplied it should be a vector giving the bounds in the following order:
#' xmin, xmax, ymin, ymax.
#'
#' @param eps A value of epsilon used in testing whether a quantity is zero,
#' mainly in the context of whether points are collinear. If anomalous errors
#' arise, it is possible that these may averted by adjusting the value of eps
#' upward or downward.
#'
#' @param normalize Should coordinates be normalized prior to calculations. If
#' \code{x} and \code{y} are in wildly different ranges it can lead to
#' tesselation and triangulation that seems off when plotted without
#' \code{\link[ggplot2]{coord_fixed}()}. Normalization of coordinates solves this.
#' The coordinates are transformed back after calculations.
#'
#' @author Thomas Lin Pedersen
#'
#' @name geom_voronoi
#' @rdname geom_voronoi
#'
#' @examples
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
	position = 'identity', na.rm = FALSE, bound = NULL, eps = 1e-9, normalize = FALSE,
	expand = 0, radius = 0, show.legend = NA, inherit.aes = TRUE, ...
) layer(
	data = data, mapping = mapping, stat = stat, geom = GeomPolygon,
	position = position, show.legend = show.legend, inherit.aes = inherit.aes,
	params = list(bound = bound, eps = eps, normalize = normalize, na.rm = na.rm, ...)
)


#' @rdname geom_voronoi
#' @importFrom scales rescale
#' @importFrom ggplot2 ggproto Stat
#' @export
StatVoronoi <- ggproto('StatVoronoi', Stat,
	compute_panel = function(data, scales, bound = NULL, eps = 1e-9, normalize = FALSE, crop = FALSE) {
		data$group <- seq_len(nrow(data))
		if (any(duplicated(data[, c('x', 'y')]))) {
			warning('stat_voronoi: dropping duplicated points', call. = FALSE)
		}
		if (normalize) {
			x_range <- range(data$x, na.rm = TRUE, finite = TRUE)
			y_range <- range(data$y, na.rm = TRUE, finite = TRUE)
			data$x <- rescale(data$x, from = x_range)
			data$y <- rescale(data$y, from = y_range)
			if (!is.null(bound)) {
				bound[1:2] <- rescale(bound[1:2], from = x_range)
				bound[3:4] <- rescale(bound[3:4], from = x_range)
			}
		}
		tiles <- to_tile(data)
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


# HELPERS -----------------------------------------------------------------


#' @importFrom sf st_multipoint st_sfc st_cast st_convex_hull st_collection_extract st_voronoi st_point st_intersects st_intersection
#' @importFrom dplyr select bind_rows %>%
to_tile <- function(data) {
	points <- data %>% select('x', 'y') %>% as.matrix() %>% st_multipoint()
	hull <- st_convex_hull(points)
	# st_voronoi returns a GEOMETRYCOLLECTION containing only polygons,
	# because a MULTIPOLYGON cannot have shared corner points.
	polys <- st_collection_extract(st_voronoi(points), 'POLYGON')
	ord <- points %>% st_sfc() %>% st_cast('POINT') %>% st_intersects(polys) %>% unlist()
	st_intersection(polys[ord], hull) %>%
		lapply(function(poly) poly %>% as.matrix() %>% as.data.frame() %>% setNames(c('x', 'y'))) %>%
		bind_rows(.id = 'group')
}

# ggplot(dm_scial, aes(DC1, DC2, fill = seq_len(ncol(scial)))) + geom_voronoi() + geom_point(shape = 21, colour = 'transparent') + scale_fill_cube_helix(r = 5, discrete = FALSE, name = 'Point')
