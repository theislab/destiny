censoring <- function(data, sigma, dists, censor_val = NULL, censor_range = NULL, missing_range = NULL, callback = invisible) {
	if (!is.null(censor_range))
		censor_range <- matrix(censor_range, ncol = 2)
	
	if (!is.null(missing_range))
		missing_range <- matrix(missing_range, ncol = 2)
	
	validate_censoring(data, sigma, dists, censor_val, censor_range, missing_range)
	
	data <- as.matrix(data)
	
	# dists is a sparse symmetrix matrix, which eigen cannot handle
	# so we just pass its filled triangle as sparse matrix and convert it back afterwards
	uplo <- dists@uplo
	half <- censoring_impl(data, sigma, structure(dists, class = 'dgCMatrix'), censor_val, censor_range, missing_range, callback)
	structure(half, class = 'dsCMatrix', uplo = uplo)
}

predict_censoring <- function(data, data2, censor_val = NULL, censor_range = NULL, missing_range = NULL, sigma) {
	if (is.null(censor_val   )) censor_val    <- NA  # this works since this will be NaN in C++ and comparison with NaN is always false
	if (is.null(censor_range )) censor_range  <- c(NA, NA)
	if (is.null(missing_range)) missing_range <- c(NA, NA)
	
	predict_censoring_impl(data, data2, censor_val, censor_range, missing_range, sigma)
}

#' @importFrom Matrix sparseMatrix
validate_censoring <- function(data, sigma, dists, censor_val, censor_range, missing_range) {
	G <- ncol(data)
	n <- nrow(data)
	
	if (any(is.na(data)) && is.null(missing_range))
		stop('Your data contains missing values (NA). You have to provide a the `missing_range` parameter.')
	
	if (is.null(censor_val) != is.null(censor_range))
		stop('You have to provide both a censoring value and a censor_range or none')
	
	if (!is.null(censor_range)) {
		if (!is.numeric(censor_val) || (length(censor_val) != 1L && length(censor_val) != G))
		stop('censor_val has to be a single numeric value, or length(censor_val) == ncol(data) must be TRUE')
		
		if (!is.numeric(censor_range) || !(nrow(censor_range) %in% c(G, 1L)) || ncol(censor_range) != 2L || any(diff(t(censor_range)) <= 0L))
		stop('censor_range has to be a numeric vector of length 2, the second of which being larger,
				 or a matrix with nrow(censor_range) == ncol(data) where each row is such a vector')
	}
	
	if (!is.null(missing_range)) {
		if (!is.numeric(missing_range) || !(nrow(missing_range) %in% c(G, 1L)) || ncol(missing_range) != 2L || any(diff(t(missing_range)) <= 0L))
		stop('missing_range has to be a numeric vector of length 2, the second of which being larger,
				 or a matrix with nrow(missing_range) == ncol(data) where each row is such a vector')
	}
	
	if (!is.numeric(sigma) || !length(sigma) %in% c(n, 1L))
		stop('sigma has to be a single numeric value or of length nrow(data)')
	
	if (!is(dists, 'dsCMatrix'))
		stop('dists has to be a dsCMatrix, not ', class(dists))
}

test_censoring <- function(censor_val, censor_range, data, missing_range) {
	if (is.null(data)) return(FALSE)
	
	has_censor_range <- !(missing(censor_range) || is.null(censor_range))
	has_missing_range <- !(missing(missing_range) || is.null(missing_range))
	if (has_censor_range || has_missing_range) return(TRUE)
	
	any(is.na(data)) || any(data == censor_val)
}
