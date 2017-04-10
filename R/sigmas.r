#' @include s4-unions.r dist-matrix-coerce.r accessor-generics.r
NULL

#' Sigmas Object
#' 
#' Holds the information about how the \code{sigma} parameter for a \link{DiffusionMap} was obtained,
#' and in this way provides a plotting function for the \link{find_sigmas} heuristic.
#' You should not need to create a Sigmas object yourself. Provide \code{sigma} to \link{DiffusionMap} instead or use \link{find_sigmas}.
#' 
#' A Sigmas object is either created by \link{find_sigmas} or by specifying the \code{sigma} parameter to \link{DiffusionMap}.
#' 
#' In the second case, if the \code{sigma} parameter is just a number,
#' the resulting \code{Sigmas} object has all slots except of \code{optimal_sigma} set to \code{NULL}.
#' 
#' @usage Sigmas(...)
#' 
#' @param ... See \dQuote{\strong{Slots}} below
#' @param object,x  \link{Sigmas} object
#' 
#' @return \code{Sigmas} creates an object of the same class
#' 
#' @slot log_sigmas     Vector of length \eqn{m} containing the \eqn{\log_{10}} of the \eqn{\sigma}s
#' @slot dim_norms      Vector of length \eqn{m-1} containing the average dimensionality \eqn{\langle p \rangle} for the respective kernel widths
#' @slot optimal_sigma  Multiple local sigmas or the mean of the two global \eqn{\sigma}s around the highest \eqn{\langle p \rangle} (\code{c(optimal_idx, optimal_idx+1L)})
#' @slot optimal_idx    The index of the highest \eqn{\langle p \rangle}.
#' @slot avrd_norms     Vector of length \eqn{m} containing the average dimensionality for the corresponding sigma.
#' 
#' @examples
#' data(guo)
#' sigs <- find_sigmas(guo, verbose = FALSE)
#' optimal_sigma(sigs)
#' print(sigs)
#' 
#' @seealso \code{\link{find_sigmas}}, the function to determine a locally optimal sigma and returning this class
#' 
#' @aliases Sigmas Sigmas-class Sigmas-methods print,Sigmas-method show,Sigmas-method optimal_sigma,Sigmas-method
#' 
#' @importFrom methods setClass
#' @name Sigmas class
#' @export Sigmas
#' @exportClass Sigmas
Sigmas <- setClass('Sigmas', slots = c(
	log_sigmas    = 'numericOrNULL',
	dim_norms     = 'numericOrNULL',
	optimal_sigma = 'numericOrNULL',
	optimal_idx   = 'integerOrNULL',
	avrd_norms    = 'numericOrNULL'))

#' @return \code{optimal_sigma} retrieves the numeric value of the optimal sigma or local sigmas
#' @name Sigmas class
#' @export
setMethod('optimal_sigma', 'Sigmas', function(object) object@optimal_sigma)

#' @name Sigmas class
#' @export
setMethod('print', 'Sigmas', function(x) {
	cat(sprintf('Sigmas (%s Steps performed)\noptimal_sigma: ', length(x@log_sigmas)))
	str(optimal_sigma(x))
	invisible(x)
})

#' @importFrom graphics plot
#' @name Sigmas class
#' @export
setMethod('show', 'Sigmas', function(object) {
	plot(object)
	invisible()
})

#' Calculate the average dimensionality for m different gaussian kernel widths (\eqn{\sigma}).
#' 
#' The sigma with the maximum value in average dimensionality is close to the ideal one.
#' Increasing step number gets this nearer to the ideal one.
#' 
#' @param data           Data set with \eqn{n} samples. Can be a \link[base]{data.frame}, \link[base]{matrix} or \link[Biobase]{ExpressionSet}.
#' @param step_size      Size of log-sigma steps
#' @param steps          Number of steps/calculations
#' @param start          Initial value to search from. (Optional. default: \eqn{\log_{10}(min(dist(data)))})
#' @param sample_rows    Number of random rows to use for sigma estimation or vector of row indices/names to use.
#'                       In the first case, only used if actually smaller than the number of available rows (Optional. default: 500)
#' @param early_exit     logical. If TRUE, return if the first local maximum is found, else keep running
#' @param ...            All parameter after this are optional and have to be specified by name
#' @param censor_val     Value regarded as uncertain. Either a single value or one for every dimension
#' @param censor_range   Uncertainity range for censoring. A length-2-vector of certainty range start and end. TODO: also allow \eqn{2\times G} matrix
#' @param missing_range  Whole data range for missing value model. Has to be specified if NAs are in the data
#' @param vars           Variables (columns) of the data to use. Specifying TRUE will select all columns (default: All floating point value columns)
#' @param verbose        logical. If TRUE, show a progress bar and plot the output
#' 
#' @return Object of class \link{Sigmas}
#' 
#' @seealso \code{\link{Sigmas}}, the class returned by this; \code{\link{DiffusionMap}}, the class this is used for
#' 
#' @examples
#' data(guo)
#' sigs <- find_sigmas(guo, verbose = TRUE)
#' DiffusionMap(guo, sigs)
#' 
#' @importFrom methods as
#' @importFrom graphics plot
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix colSums rowSums
#' @importFrom proxy dist
#' @export
find_sigmas <- function(
	data,
	step_size = .1, steps = 10L,
	start = NULL,
	sample_rows = 500L,
	early_exit = FALSE,
	...,
	censor_val = NULL, censor_range = NULL,
	missing_range = NULL,
	vars = NULL,
	verbose = TRUE
) {
	data <- extract_doublematrix(data, vars)
	
	if (any(is.na(data)))
		data <- as.matrix(hotdeck(data, imp_var = FALSE))
	
	stopifnot(steps >= 3L)
	
	if (length(sample_rows) > 1L) {
		data <- data[sample_rows, ]
	} else if (nrow(data) > sample_rows) {
		sample_idx <- sample(nrow(data), sample_rows)
		data <- data[sample_idx, ]
	}
	
	n <- nrow(data)
	
	dists <- dist(data)
	
	min_dist <- min(dists)
	if (min_dist == 0)
		stop('Minimum distance in the data may not be 0')
	dists <- as(dists, 'symmetricMatrix')
	
	if (is.null(start))
		start <- log10(min_dist)
	#if (missing(step_size))
	#  step_size = min_dist / steps
	
	if (verbose) print(c(min_dist = min_dist, start = start, step_size = step_size))
	
	get_trans_p <-
		if (test_censoring(censor_val, censor_range, data, missing_range)) {
			function(sigma) censoring(data, sigma, dists, censor_val, censor_range, missing_range)
		} else {
			msqd <- -(dists ^ 2)
			function(sigma) exp(msqd / (2 * sigma ^ 2))
		}
	
	do_step <- function(i) {
		# i can be negative!
		log_sigma <- start + i*step_size
		trans_p <- get_trans_p(10 ^ log_sigma)
		
		diag_d <- colSums(trans_p, na.rm = TRUE)
		
		list(avrd_norm = (sum(log10(diag_d/n) / diag_d)) / sum(1 / diag_d),
		     log_sigma = log_sigma)
	}
	
	avrd_norms <- numeric(steps)
	log_sigmas <- numeric(steps)
	dim_norms  <- numeric(steps - 1)
	
	step_diff <- function(step) {
		idxs <- c(step, step - 1)
		diff(avrd_norms[idxs]) / diff(log_sigmas[idxs])
	}
	
	a0 <- do_step(0L)
	avrd_norms[[1L]] <- a0$avrd_norm
	log_sigmas[[1L]] <- a0$log_sigma
	
	a1 <- do_step(1)
	dir <- 1L
	avrd_norms[[2L]] <- a1$avrd_norm
	log_sigmas[[2L]] <- a1$log_sigma
	
	if (step_diff(2L) < 0) {
		a1 <- do_step(-1L)
		dir <- -1L
		avrd_norms[[2L]] <- a1$avrd_norm
		log_sigmas[[2L]] <- a1$log_sigma
	}
	
	dim_norms[[1L]] <- step_diff(2L)
	
	if (verbose) pb <- txtProgressBar(2L, steps, 1L, style = 3)
	for (step in seq(2L, steps)) {
		a_i = do_step(dir * (step - 1L))
		avrd_norms[[step]] <- a_i$avrd_norm
		log_sigmas[[step]] <- a_i$log_sigma
		
		dif_step <- step - 1
		dim_norms[[dif_step]] <- step_diff(step)
		
		if (verbose) setTxtProgressBar(pb, step)
		
		if (early_exit && step > 2 && dim_norms[[dif_step]] < dim_norms[[dif_step - 1L]]) {
			avrd_norms <- avrd_norms[seq_len(step)]
			log_sigmas <- log_sigmas[seq_len(step)]
			dim_norms  <- dim_norms[seq_len(dif_step)]
			break
		}
	}
	if (verbose) {
		setTxtProgressBar(pb, steps)
		close(pb)
	}
	
	if (early_exit && step == steps) warning('All steps were exhausted without finding a maximum. Using last encountered sigma')
	
	optimal_idx <- which.max(dim_norms)
	
	ret <- Sigmas(
		log_sigmas = log_sigmas,
		dim_norms  = dim_norms,
		optimal_sigma = 10 ^ mean(log_sigmas[c(optimal_idx, optimal_idx + 1L)]),
		optimal_idx   = optimal_idx,
		avrd_norms = avrd_norms)
	
	if (verbose) plot(ret)
	
	ret
}
