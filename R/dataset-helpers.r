#' @importFrom Biobase exprs
extract_doublematrix <- function(data, vars = NULL) {
	if (is.data.frame(data)) {
		data <- as.matrix(data[, sapply(data, is.double)])
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


#' @importFrom Biobase sampleNames
n_samples <- function(data, distances) {
	if (is.null(data)) nrow(distances)
	else if (is(data, 'ExpressionSet')) length(sampleNames(data))
	else nrow(data)
}


#' @importFrom Biobase featureNames
n_features <- function(data, distances = NULL) {
	if (is.null(data)) ncol(distances)
	else if (is(data, 'ExpressionSet')) length(featureNames(data))
	else ncol(data)
}
