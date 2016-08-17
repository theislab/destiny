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
