#' @include dpt.r
NULL

#' Plot DPT
#' 
#' Plots diffusion components from a Diffusion Map and the accompanying Diffusion Pseudo Time (\code{\link{DPT}})
#' 
#' @param x           A \code{\link{DPT}} object.
#' @param y,root      Root branch ID. Will be used as the start of the DPT.
#'                    (If longer than size 1, will be interpreted as \code{c(root, branches)})
#' @param branches    Numeric Branch IDs. Are used as target(s) for the path(s) to draw.
#' @param dcx,dcy     The dimensions to use from the DiffusionMap
#' @param divide      If \code{col_by = 'branch'}, this specifies which branches to divide. (see \code{\link{branch_divide}})
#' @param w_width     Window width for smoothing the path (see \code{\link[smoother]{smth.gaussian}})
#' @param col_by      Color by 'dpt' (DPT starting at \code{branches[[1]]}), 'branch', or a veriable of the data.
#' @param col_path    Colors for the path or a function creating n colors
#' @param col_tip     Color for branch tips
#' @param ...         Graphical parameters supplied to \code{\link{plot.DiffusionMap}}
#' @param pal         Palette to use for coloring the points
#' 
#' @aliases plot,DPT,numeric-method plot,DPT,missing-method
#' 
#' @examples
#' data(guo_norm)
#' dm <- DiffusionMap(guo_norm)
#' dpt <- DPT(dm, branching = TRUE)
#' plot(dpt)
#' plot(dpt, 2L,      col_by = 'branch')
#' plot(dpt, 1L, 2:3, col_by = 'num_cells')
#' 
#' @importFrom graphics plot points
#' @importFrom methods is setMethod
#' @importFrom scales colour_ramp rescale
#' @export
plot.DPT <- function(
	x, root = 1L,
	branches = integer(0L),
	dcs = 1:2,
	divide = integer(0L),
	w_width = .1,
	col_by = 'dpt',
	col_path = palette(),
	col_tip = 'red',
	...,
	pal = switch(col_by, dpt = cube_helix, palette())
) {
	dpt <- x
	root     <- as.integer(root)
	branches <- as.integer(branches)
	if (length(root) < 1L) stop('root needs to be specified')
	if (length(root) > 1L && length(branches) > 0L)
		stop('(length(root), length(branches)) needs to be (1, 0-n) or (2-n, 0), but is (', length(root), ', ', length(branches), ')')
	stopifnot(length(dcs) %in% 2:3)
	
	if (length(root) > 1L && length(branches) == 0L) {
		branches <- root[-1]
		root <- root[[1]]
	}
	
	evs <- eigenvectors(dpt@dm)[, dcs]
	pt_vec <- dpt@dpt[, root]
	
	dpt_flat <- branch_divide(dpt, divide)
	
	plot_more <- function(p) {
		for (b in seq_along(branches)) {
			idx <- dpt@branch[, 1] %in% c(root, branches[[b]])
			path <- average_path(pt_vec[idx], evs[idx, ], w_width)
			points(path[, 1L], path[, 2L], 'l', col = col_path[[b]])
		}
		
		tips <- evs[dpt_flat@tips[, 1], ]
		if (is.null(p)) {  # 2d plot
			points(tips, col = col_tip)
		} else if (is.list(p)) {  # scatterplot3d
			p$points3d(tips, col = col_tip)
		} else if (is.vector(p)) {  # rgl
			rgl::points3d(tips, col = col_tip)
		} else stop('unknown p passed to plot_more (class: ', class(p), ')')
	}
	
	args <- switch(col_by,
		dpt    = list(col = pt_vec,               draw_legend = TRUE, legend_main = 'DPT'),
		branch = list(col = dpt_flat@branch[, 1], draw_legend = TRUE, legend_main = 'Branch'),
		         list(col_by = col_by))
	
	do.call(plot, c(list(dpt@dm, dcs, plot_more = plot_more, pal = pal), args, list(...)))
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
