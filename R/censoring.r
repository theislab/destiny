censoring <- function(data, censor_val = NULL, censor_range = NULL, missing_range = NULL, sigma, nns = NULL, callback = function(i) {}) {
	if (!is.null(censor_range))
		censor_range <- matrix(censor_range, ncol = 2)
	
	if (!is.null(missing_range))
		missing_range <- matrix(missing_range, ncol = 2)
	
	validate_censoring(censor_val, censor_range, missing_range, data, sigma, nns)
	
	data <- as.matrix(data)
	
	censoring_impl(data, censor_val, censor_range, missing_range, sigma, nns, callback)
}

predict_censoring <- function(data, data2, censor_val = NULL, censor_range = NULL, missing_range = NULL, sigma) {
	if (is.null(censor_val   )) censor_val    <- NA  # this works since this will be NaN in C++ and comparison with NaN is always false
	if (is.null(censor_range )) censor_range  <- c(NA, NA)
	if (is.null(missing_range)) missing_range <- c(NA, NA)
	
	predict_censoring_impl(data, data2, censor_val, censor_range, missing_range, sigma)
}

validate_censoring <- function(censor_val, censor_range, missing_range, data, sigma, nns) {
	G <- ncol(data)
	no_censor_val     <- missing(censor_val)    || is.null(censor_val)
	no_censor_range   <- missing(censor_range)  || is.null(censor_range)
	no_missing_range  <- missing(missing_range) || is.null(missing_range)
	
	if (any(is.na(data)) && no_missing_range)
		stop('Your data contains missing values (NA). You have to provide a the `missing_range` parameter.')
	
	if (no_censor_val != no_censor_range)
		stop('You have to provide both a censoring value and a censor_range or none')
	
	if (!no_censor_range) {
		if (!is.numeric(censor_val) || (length(censor_val) != 1 && length(censor_val) != G))
		stop('censor_val has to be a single numeric value, or length(censor_val) == ncol(data) must be TRUE')
		
		if (!is.numeric(censor_range) || !(nrow(censor_range) %in% c(G, 1)) || ncol(censor_range) != 2 || any(diff(t(censor_range)) <= 0))
		stop('censor_range has to be a numeric vector of length 2, the second of which being larger,
				 or a matrix with nrow(censor_range) == ncol(data) where each row is such a vector')
	}
	
	if (!no_missing_range) {
		if (!is.numeric(missing_range) || !(nrow(missing_range) %in% c(G, 1)) || ncol(missing_range) != 2 || any(diff(t(missing_range)) <= 0))
		stop('missing_range has to be a numeric vector of length 2, the second of which being larger,
				 or a matrix with nrow(missing_range) == ncol(data) where each row is such a vector')
	}
	
	if (!is.numeric(sigma) || length(sigma) != 1)
		stop('sigma has to be a single numeric value (local sigma with censoring not yet supported)')
	
	if (!missing(nns) && !is.null(nns) && !is.integer(nns))
		stop('nns has to be a integer matrix or NULL')
}

test_censoring <- function(censor_val, censor_range, data, missing_range) {
	!(missing(censor_range) || is.null(censor_range)) ||
		any(is.na(data)) || any(data == censor_val) ||
		!(missing(missing_range) || is.null(missing_range))
}
