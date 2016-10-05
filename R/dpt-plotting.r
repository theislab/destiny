#' @include dpt.r
NULL

#' Plot DPT
#' 
#' Plots diffusion components from a Diffusion Map and the accompanying Diffusion Pseudo Time (\code{\link{DPT}})
#' 
#' @param x           A \code{\link{DPT}} object.
#' @param y,root      Root branch ID. Will be used as the start of the DPT.
#'                    (If longer than size 1, will be interpreted as \code{c(root, branches)})
#' @param paths_to    Numeric Branch IDs. Are used as target(s) for the path(s) to draw.
#' @param dcs         The dimensions to use from the DiffusionMap
#' @param divide      If \code{col_by = 'branch'}, this specifies which branches to divide. (see \code{\link{branch_divide}})
#' @param w_width     Window width for smoothing the path (see \code{\link[smoother]{smth.gaussian}})
#' @param col_by      Color by 'dpt' (DPT starting at \code{branches[[1]]}), 'branch', or a veriable of the data.
#' @param col_path    Colors for the path or a function creating n colors
#' @param col_tip     Color for branch tips
#' @param ...         Graphical parameters supplied to \code{\link{plot.DiffusionMap}}
#' @param col         See \code{\link{plot.DiffusionMap}}. This overrides \code{col_by}
#' @param legend_main See \code{\link{plot.DiffusionMap}}.
#' 
#' @return The return value of the underlying call is returned, i.e. a scatterplot3d or rgl object for 3D plots.
#' 
#' @aliases plot.DPT plot,DPT,numeric-method plot,DPT,missing-method
#' 
#' @examples
#' data(guo_norm)
#' dm <- DiffusionMap(guo_norm)
#' dpt <- DPT(dm, branching = TRUE)
#' plot(dpt)
#' plot(dpt, 2L,      col_by = 'branch')
#' plot(dpt, 1L, 2:3, col_by = 'num_cells')
#' plot(dpt, col_by = 'DPT3')
#' 
#' @importFrom graphics plot points
#' @importFrom methods is setMethod
#' @importFrom scales colour_ramp rescale
#' @importFrom utils capture.output
#' @export
plot.DPT <- function(
	x, root = 1L,
	paths_to = integer(0L),
	dcs = 1:2,
	divide = integer(0L),
	w_width = .1,
	col_by = 'dpt',
	col_path = rev(palette()),
	col_tip = 'red',
	...,
	col = NULL,
	legend_main = col_by
) {
	dpt      <- x
	root     <- as.integer(root)
	paths_to <- as.integer(paths_to)
	if (length(root) < 1L) stop('root needs to be specified')
	if (length(root) > 1L && length(paths_to) > 0L)
		stop('(length(root), length(paths_to)) needs to be (1, 0-n) or (2-n, 0), but is (', length(root), ', ', length(paths_to), ')')
	stopifnot(length(dcs) %in% 2:3)
	
	if (length(root) > 1L && length(paths_to) == 0L) {
		paths_to <- root[-1]
		root <- root[[1]]
	}
	
	evs <- flipped_dcs(dpt@dm, dcs)
	
	dpt_flat <- branch_divide(dpt, divide)
	
	pt_vec <- dpt_for_branch(dpt_flat, root)
	
	plot_paths <- function(p, ..., rescale) {
		plot_points <- get_plot_fn(p)
		rescale_fun <-
			if (is.null(rescale)) identity
			else function(x) rescale_mat(x, from = rescale$from, to = rescale$to)
		
		for (b in seq_along(paths_to)) {
			idx <- dpt@branch[, 1] %in% c(root, paths_to[[b]])
			path <- average_path(pt_vec[idx], evs[idx, ], w_width)
			plot_points(rescale_fun(path), type = 'l', col = col_path[[b]], ...)
		}
		
		tips <- evs[dpt_flat@tips[, 1], ]
		plot_points(rescale_fun(tips), col = col_tip, ...)
	}
	
	col <-
		if (!is.null(col)) col
		else switch(col_by,
			dpt    = pt_vec,
			branch = ,
			Branch = dpt_flat@branch[, 1],
			dpt[[col_by]])
	
	legend_main <- switch(legend_main, dpt = 'DPT', branch = 'Branch', legend_main)
	
	args <- list(
		dpt@dm, dcs,
		plot_more = plot_paths,
		legend_main = legend_main,
		col = col,
		...)
	
	if (!identical(Sys.getenv('LOG_LEVEL'), '')) message('Args:\n', paste(capture.output(print(args)), collapse = '\n'))
	do.call(plot, args)
}

#' @name plot.DPT
#' @export
setMethod('plot', c('DPT', 'numeric'), function(x, y, ...) plot.DPT(x, y, ...))

#' @name plot.DPT
#' @export
setMethod('plot', c('DPT', 'missing'), function(x, y, ...) {
	args <- list(...)
	root <- if (is.null(args$root)) 1L else args$root
	args$root <- NULL
	
	do.call(plot.DPT, c(list(x, root), args))
})


#' @importFrom graphics plot
#' @importFrom smoother smth.gaussian
average_path <- function(pt, x, w_width = .1) {
	stopifnot(identical(nrow(x), length(pt)))
	as.data.frame(apply(x[order(pt), ], 2, function(col) smth.gaussian(col, w_width, tails = TRUE)))
}


get_plot_fn <- function(p) {
	if (is.null(p)) {  # 2d plot
		points
	} else if (is.list(p) && 'points3d' %in% names(p)) {# scatterplot3d
		p$points3d
	} else if (is(p, 'rglHighlevel')) {  # rgl
		function(x, y = NULL, z = NULL, type = 'p', ...)
			switch(type, p = rgl::points3d, l = rgl::lines3d, stop)(x, y, z, ...)
	} else stop('unknown p passed to plot_more (class: ', class(p), ')')
}
