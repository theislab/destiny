censoring <- function(data, censor.val = NULL, censor.range = NULL, missing.range = NULL, sigma, nns = NULL, callback = function(i) {}) {
	if (!is.null(censor.range))
		censor.range <- matrix(censor.range, ncol = 2)
	
	if (!is.null(missing.range))
		missing.range <- matrix(missing.range, ncol = 2)
	
	validate.censoring(censor.val, censor.range, missing.range, data, sigma, nns)
	
	data <- as.matrix(data)
	
	censoring_impl(data, censor.val, censor.range, missing.range, sigma, nns, callback)
}

predict.censoring <- function(data, data2, censor.val = NULL, censor.range = NULL, missing.range = NULL, sigma) {
	if (is.null(censor.val   )) censor.val    <- NA  # this works since this will be NaN in C++ and comparison with NaN is always false
	if (is.null(censor.range )) censor.range  <- c(NA, NA)
	if (is.null(missing.range)) missing.range <- c(NA, NA)
	
	predict_censoring_impl(data, data2, censor.val, censor.range, missing.range, sigma)
}

validate.censoring <- function(censor.val, censor.range, missing.range, data, sigma, nns) {
	G <- ncol(data)
	no.censor.val     <- missing(censor.val)    || is.null(censor.val)
	no.censor.range   <- missing(censor.range)  || is.null(censor.range)
	no.missing.range  <- missing(missing.range) || is.null(missing.range)
	
	if (any(is.na(data)) && no.missing.range)
		stop('Your data contains missing values (NA). You have to provide a the `missing.range` parameter.')
	
	if (no.censor.val != no.censor.range)
		stop('You have to provide both a censoring value and a censor.range or none')
	
	if (!no.censor.range) {
		if (!is.numeric(censor.val) || (length(censor.val) != 1 && length(censor.val) != G))
		stop('censor.val has to be a single numeric value, or length(censor.val) == ncol(data) must be TRUE')
		
		if (!is.numeric(censor.range) || !(nrow(censor.range) %in% c(G, 1)) || ncol(censor.range) != 2 || any(diff(t(censor.range)) <= 0))
		stop('censor.range has to be a numeric vector of length 2, the second of which being larger,
				 or a matrix with nrow(censor.range) == ncol(data) where each row is such a vector')
	}
	
	if (!no.missing.range) {
		if (!is.numeric(missing.range) || !(nrow(missing.range) %in% c(G, 1)) || ncol(missing.range) != 2 || any(diff(t(missing.range)) <= 0))
		stop('missing.range has to be a numeric vector of length 2, the second of which being larger,
				 or a matrix with nrow(missing.range) == ncol(data) where each row is such a vector')
	}
	
	if (!is.numeric(sigma) || length(sigma) != 1)
		stop('sigma has to be a single numeric value')
	
	if (!missing(nns) && !is.null(nns) && !is.integer(nns))
		stop('nns has to be a integer matrix or NULL')
}

test.censoring <- function(censor.val, censor.range, data, missing.range) {
	!(missing(censor.range) || is.null(censor.range)) ||
		any(is.na(data)) || any(data == censor.val) ||
		!(missing(missing.range) || is.null(missing.range))
}

no_censoring_slow <- function(sigma, nn.index, nn.dist, cb) {
	k <- ncol(nn.index)
	n <- nrow(nn.index)
	stopifnot(k == ncol(nn.dist))
	stopifnot(n == nrow(nn.dist))
	
	trans.p <- sparseMatrix(NULL, NULL, x = numeric(0), dims = c(n, n))
	
	for (i in seq_len(k)) {
		pairs.i <- cbind(seq_len(n), nn.index[, i])
		trans.p[pairs.i] <- exp(-nn.dist[, i] ^ 2 / (2 * sigma ^ 2))
		cb(i)
	}
	
	trans.p
}
