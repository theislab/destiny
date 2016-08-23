#' @include dpt.r
NULL

#' Plot DPT
#' 
#' Plots diffusion components from a Diffusion Map and the accompanying Diffusion Pseudo Time (\code{\link{DPT}})
#' 
#' @param x           A \code{\link{DPT}} object.
#' @param y,branches  Numeric Branch IDs. The first one will be used as the start of the DPT, subsequent ones for the path(s).
#' @param dcx,dcy     The dimensions to use from the DiffusionMap
#' @param w_width     Window width for smoothing the path (see \code{\link[smoother]{smth.gaussian}})
#' @param col_by      Color by 'dpt' (DPT starting at \code{branches[[1]]}), 'branch', or a veriable of the data.
#' @param col_map     Color Map for DPT coloring. Sequence of colors or a function accepting the number of colors to create
#' @param col_path    Colors for the path or a function creating n colors
#' @param col_tip     Color for branch tips
#' @param ...         Graphical parameters supplied to \code{\link{plot.DiffusionMap}}
#' 
#' @aliases plot,DPT,numeric-method plot,DPT,missing-method
#' 
#' @examples
#' data(guo_norm)
#' dm <- DiffusionMap(guo_norm)
#' dpt <- DPT(dm, branching = TRUE)
#' plot(dpt)
#' plot(dpt, 2L,  col_by = 'branch')
#' plot(dpt, 1:3, col_by = 'num_cells')
#' 
#' @importFrom graphics plot points
#' @importFrom methods is setMethod
#' @importFrom scales colour_ramp rescale
#' @export
plot.DPT <- function(
	x, branches = 1L,
	dcx = 1L, dcy = 2L,
	w_width = .1,
	col_by = 'dpt',
	col_map = cube_helix,
	col_path = palette(),
	col_tip = 'red',
	...
) {
	dpt <- x
	branches <- as.integer(branches)
	stopifnot(length(dcx) == 1L, length(dcy) == 1L)
	
	if (!is.function(col_map)) {
		col_map <- colour_ramp(col_map)
	}
	
	evs <- eigenvectors(dpt@dm)
	pt_vec <- dpt@dpt[, branches[[1L]]]
	
	plot_more <- function() {
		for (b in seq_along(branches[-1])) {
			idx <- dpt@branch[, 1] %in% c(branches[[1]], branches[[b + 1L]])
			path <- average_path(pt_vec[idx], evs[idx, c(dcx, dcy)], w_width)
			points(path[, 1L], path[, 2L], 'l', col = col_path[[b]])
		}
		
		points(evs[dpt@tips[, 1], dcx], evs[dpt@tips[, 1], dcy], col = col_tip)
	}
	
	args <- switch(col_by,
		dpt    = list(col = pt_vec, pal = col_map, draw_legend = TRUE, legend_main = 'DPT'),
		branch = list(col = dpt@branch[, 1],       draw_legend = TRUE, legend_main = 'Branch'),
		         list(col_by = col_by))
	
	do.call(plot, c(list(dpt@dm, c(dcx, dcy), plot_more = plot_more), args, ...))
}

#' @name plot.DPT
#' @export
setMethod('plot', c('DPT', 'numeric'), function(x, y, ...) plot.DPT(x, y, ...))

#' @name plot.DPT
#' @export
setMethod('plot', c('DPT', 'missing'), function(x, y, ...) plot.DPT(x, 1L, ...))


#' @importFrom graphics plot
#' @importFrom smoother smth.gaussian
average_path <- function(pt, x, w_width = .1) {
	stopifnot(identical(nrow(x), length(pt)))
	as.data.frame(apply(x[order(pt), ], 2, function(col) smth.gaussian(col, w_width, tails = TRUE)))
}
