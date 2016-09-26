#' @importFrom Biobase exprs
extract_doublematrix <- function(data, vars = NULL) {
	if (is.data.frame(data)) {
		data <- as.matrix(data[sapply(data, is.double)])
	} else if (inherits(data, 'ExpressionSet')) {
		data <- t(exprs(data))
	} else if (!is.matrix(data)) {
		stop('Data needs to be matrix, data.frame or ExpressionSet')
	}
	dupes <- duplicated(data)
	if (any(dupes)) {
		data <- data[!dupes, ]
		warning('Duplicate rows removed from data. Consider explicitly using `df[!duplicated(df), ]`')
	}
	
	if (!is.null(vars))
		data <- data[, vars]
	data
}


stopifsmall <- function(max_dist) {
	if (max_dist < .Machine$double.eps)
		stop(sprintf(
			'The supplied sigma is not large enough. Please select a larger one.
			find_sigmas(data) should return one with the right order of magnitude. (max dist. is %.3e)',
			max_dist))
}


verbose_timing <- function(verbose, msg, expr) {
	if (verbose) {
		cat(sprintf('%s...', msg))
		dif <- system.time({
			r <- force(expr)
		})
		cat(sprintf('...done. Time: %.2fs\n', dif[['elapsed']]))
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


rescale_mat <- function(mat, ...) apply(mat, 2L, scales::rescale, ...)
