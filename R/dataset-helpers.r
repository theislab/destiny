#' @importFrom Biobase exprs
dataset_extract_doublematrix <- function(data, vars = NULL) {
	if (is.data.frame(data)) {
		if (is.null(vars)) {
			data <- data[, sapply(data, is.double)]
		} else {  # If we have a data frame subset early, before matrix conversion
			data <- data[, vars]
			vars <- NULL
		}
		data <- as.matrix(data)
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
dataset_n_observations <- function(data, distances) {
	if (is.null(data)) nrow(distances)
	else if (is(data, 'ExpressionSet')) length(sampleNames(data))
	else nrow(data)
}


#' @importFrom Biobase featureNames
dataset_n_features <- function(data, distances = NULL, vars = NULL) {
	if (is.null(data)) ncol(distances)
	else if (is(data, 'ExpressionSet')) length(featureNames(data))
	else if (!is.null(vars)) ncol(data[, vars])
	else if (is.data.frame(data)) ncol(data[, sapply(data, is.double)])
	else ncol(data)
}


#' @importFrom methods canCoerce
dataset_to_df <- function(dta, row.names = NULL, optional = FALSE, ...) {
	# The ExpressionSet as.data.frame sucks
	if (is(dta, 'ExpressionSet')) {
		cbind(as.data.frame(t(exprs(dta)), row.names, optional, ...), pData(dta))
	} else if (canCoerce(dta, 'data.frame')) {
		as(dta, 'data.frame')
	} else if (!is.null(getS3method('as.data.frame', class(data)[[1L]], optional = TRUE))) {
		as.data.frame(dta, row.names, optional, ...)
	} else NULL
}


dataset_names <- function(dta) {
	if (is(dta, 'ExpressionSet')) {
		c(featureNames(dta), varLabels(dta))
	} else {
		colnames(dta)
	}
}


#' @importFrom Biobase featureNames exprs
dataset_get_feature <- function(dta, f) {
	force(f)  # else dta[, f] would be dta[, ] and valid
	tryCatch({
		if (is(dta, 'ExpressionSet') && f %in% featureNames(dta)) {
			exprs(dta)[f, ]
		} else if (is(dta, 'ExpressionSet') || is.data.frame(dta)) {
			# data.frame or phenoData
			dta[[f]]
		} else {
			dta[, f]
		}
	}, error = function(e) stop(sprintf('Invalid `%s`: No column, annotation, or feature found with name %s', deparse(substitute(f)), dQuote(f))))
}
