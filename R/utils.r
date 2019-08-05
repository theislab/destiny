# There is currently a problem on windows which can't display chars in th
# text/plain output, which are not available in the current locale.
# See https://github.com/IRkernel/repr/issues/28#issuecomment-208574856
#' @importFrom utils capture.output
char_fallback <- function(char, default) {
	real_len <- nchar(char)
	r_len <- nchar(capture.output(cat(char)))
	if (real_len == r_len) char else default
}
sym_times <- char_fallback('\u00D7', 'x')


stopifsmall <- function(max_dist) {
	if (max_dist < .Machine$double.eps)
		stop(sprintf(
			'The supplied sigma is not large enough. Please select a larger one.
			find_sigmas(data) should return one with the right order of magnitude. (max dist. is %.3e)',
			max_dist))
}


stopifparams <- function(...) {
	args <- substitute(...())  # Unevaluated arguments
	if (length(args) == 0) return(invisible())
	nms <- if (is.null(names(args))) rep_len('', length(args)) else names(args)
	args_str <- paste0(nms, ifelse(nms == '', '', ' = '), sapply(args, deparse))
	stop('Unused argument(s) (', paste(args_str, collapse = ', '), ')', call. = FALSE)
}


#' @importFrom utils flush.console
verbose_timing <- function(verbose, msg, expr) {
	if (verbose) {
		cat(sprintf('%s...', msg))
		flush.console()
		dif <- system.time({
			r <- force(expr)
		})
		cat(sprintf('...done. Time: %.2fs\n', dif[['elapsed']]))
		flush.console()
		r
	} else expr
}


#' @importFrom Matrix Diagonal
#' @importMethodsFrom Matrix solve
accumulated_transitions <- function(dm) {
	if (!is.null(dm@data_env$propagations)) {  # compat
		dm@data_env$accumulated_transitions <- dm@data_env$propagations
		rm('propagations', envir = dm@data_env)
	}
	
	if (is.null(dm@data_env$accumulated_transitions)) {
		if (is.null(dm@transitions))
			stop('DiffusionMap was created with suppress_dpt = TRUE')
		
		n <- length(dm@d_norm)
		
		phi0 <- dm@d_norm / sqrt(sum(dm@d_norm ^ 2))
		inv <- solve(Diagonal(n) - dm@transitions + phi0 %*% t(phi0))
		dm@data_env$accumulated_transitions <- inv - Diagonal(n)
	}
	
	dm@data_env$accumulated_transitions
}


hasattr <- function(x, which) !is.null(attr(x, which, exact = TRUE))


flipped_dcs <- function(d, dcs) {
	if (is(d, 'DiffusionMap')) d <- eigenvectors(d)
	
	evs <- as.matrix(d[, abs(dcs)])
	evs[, dcs < 0] <- -evs[, dcs < 0]
	evs
}


rescale_mat <- function(mat, rescale) {
	if (is.list(rescale)) {
		stopifnot(setequal(dimnames(rescale), c('from', 'to')))
		rv <- apply(mat, 2L, scales::rescale, rescale$to, rescale$from)
	} else if (is.array(rescale)) {
		stopifnot(length(dim(rescale)) == 3L)
		stopifnot(ncol(mat) == ncol(rescale))
		stopifnot(dim(rescale)[[1L]] == 2L)
		stopifnot(dim(rescale)[[3L]] == 2L)
		
		col_type <- get(typeof(mat))
		rv <- vapply(seq_len(ncol(mat)), function(d) {
			scales::rescale(mat[, d], rescale['to', d, ], rescale['from', d, ])
		}, col_type(nrow(mat)))
	}
	stopifnot(all(dim(rv) == dim(mat)))
	dimnames(rv) <- dimnames(mat)
	rv
}


runs <- function(vec) {
	enc <- rle(vec)
	enc$values <- make.unique(enc$values, '_')
	inverse.rle(enc)
}


upper.tri.sparse <- function(x, diag = FALSE) {
	# Works just like upper.tri() but doesn't forcibly coerce large 'sparseMatrix' back to 'matrix'
	if (diag)
		row(x) <= col(x)
	else row(x) < col(x)
}


get_louvain_clusters <- function(transitions) {
	graph <- igraph::graph_from_adjacency_matrix(transitions, 'undirected', weighted = TRUE)
	as.integer(unclass(igraph::membership(igraph::cluster_louvain(graph))))
}


#' @importFrom BiocGenerics duplicated
setMethod('duplicated', 'dgCMatrix', function(x, incomparables = FALSE, MARGIN = 1L, ...) {
	MARGIN <- as.integer(MARGIN)
	n <- nrow(x)
	p <- ncol(x)
	j <- rep(seq_len(p), diff(x@p))
	i <- x@i + 1
	v <- x@x
	
	if (MARGIN == 1L) {  # rows
		names(v) <- j
		splits <- split(v, i)
		is_empty <- setdiff(seq_len(n), i)
	} else if (MARGIN == 2L) {  # columns
		names(v) <- i
		splits <- split(v, j)
		is_empty <- setdiff(seq_len(p), j)
	} else stop('Invalid MARGIN ', MARGIN, ', matrices only have rows (1) and columns (2).')
	
	result <- duplicated.default(splits)
	if (!any(is_empty)) return(result)
	
	out <- logical(if (MARGIN == 1L) n else p)
	out[-is_empty] <- result
	if (length(is_empty) > 1)
		out[is_empty[-1]] <- TRUE
	out
})

if (!exists('hcl.colors', 'package:grDevices'))
	hcl.colors <- function(
		n,
		palette = 'viridis',
		alpha = NULL,
		rev = FALSE,
		fixup = TRUE
	) {
		option <- switch(
			palette,
			magma = 'A', inferno = 'B', plasma = 'C', viridis = 'D',
			stop('You are using an old R version and the hcl.colors replacement will not work'))
		viridisLite::viridis(n, alpha, direction = if (rev) -1 else 1, option = option)
	}

