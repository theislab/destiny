#' @import Matrix
#' @include s4-null-unions.r dist-matrix-coerce.r
NULL

#' Sigmas Object
#' 
#' Holds the information about how the \code{sigma} parameter for a \link{DiffusionMap} was obtained,
#' and in this way provides a plotting function for the \link{find.sigmas} heuristic.
#' You should not need to create a Sigmas object yourself. Provide \code{sigma} to \link{DiffusionMap} instead or use \link{find.sigmas}.
#' 
#' A Sigmas object is either created by \link{find.sigmas} or by specifying the \code{sigma} parameter to \link{DiffusionMap}.
#' 
#' In the second case, if the \code{sigma} parameter is just a number,
#' the resulting \code{Sigmas} object has all slots except of \code{optimal.sigma} set to \code{NULL}.
#' 
#' @usage Sigmas(...)
#' 
#' @param ... See \dQuote{\strong{Slots}} below
#' @param object,x  \link{Sigmas} object
#' 
#' @return \code{Sigmas} creates an object of the same class
#' 
#' @slot log.sigmas     Vector of length \eqn{m} containing the \eqn{\log_{10}} of the \eqn{\sigma}s
#' @slot dim.norms      Vector of length \eqn{m-1} containing the average dimensionality \eqn{\langle p \rangle} for the respective kernel widths
#' @slot optimal.sigma  Mean of the two \eqn{\sigma}s around the highest \eqn{\langle p \rangle} (\code{c(optimal.idx, optimal.idx+1L)})
#' @slot optimal.idx    The index of the highest \eqn{\langle p \rangle}.
#' @slot avrd.norms     Vector of length \eqn{m} containing the average dimensionality for the corresponding sigma.
#' 
#' @examples
#' data(guo)
#' sigs <- find.sigmas(guo, verbose = FALSE)
#' optimal.sigma(sigs)
#' print(sigs)
#' 
#' @seealso \code{\link{find.sigmas}}, the function to determine a locally optimal sigma and returning this class
#' 
#' @aliases Sigmas Sigmas-class Sigmas-methods print,Sigmas-method show,Sigmas-method optimal.sigma optimal.sigma,Sigmas-method
#' @name Sigmas class
#' @export Sigmas
#' @exportClass Sigmas
Sigmas <- setClass('Sigmas', slots = c(
	log.sigmas    = 'numericOrNULL',
	dim.norms     = 'numericOrNULL',
	optimal.sigma = 'numeric',
	optimal.idx   = 'integerOrNULL',
	avrd.norms    = 'numericOrNULL'))

#' @return \code{optimal.sigma} retrieves the numeric value of the optimal sigma
#' 
#' @name Sigmas class
#' @export
setGeneric('optimal.sigma', function(object) standardGeneric('optimal.sigma'))

#' @name Sigmas class
#' @export
setMethod('optimal.sigma', 'Sigmas', function(object) object@optimal.sigma)

#' @name Sigmas class
#' @export
setMethod('print', 'Sigmas', function(x) {
	cat(sprintf('Sigmas (%s Steps performed)
optimal.sigma: %s',
		length(x@log.sigmas),
		optimal.sigma(x)))
	invisible(x)
})

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
#' @param step.size      Size of log-sigma steps
#' @param steps          Number of steps/calculations
#' @param start          Initial value to search from. (Optional. default: \eqn{\log_10(min(dist(data)))})
#' @param sample.rows    Number of random rows to use for sigma estimation or vector of row indices/names to use.
#'                       In the first case, only used if actually smaller than the number of available rows (Optional. default: 500)
#' @param early.exit     logical. If TRUE, return if the first local maximum is found, else keep running
#' @param ...            All parameter after this are optional and have to be specified by name
#' @param censor.val     Value regarded as uncertain. Either a single value or one for every dimension
#' @param censor.range   Uncertainity range for censoring. A length-2-vector of certainty range start and end. TODO: also allow \eqn{2\times G} matrix
#' @param missing.range  Whole data range for missing value model. Has to be specified if NAs are in the data
#' @param vars           Variables (columns) of the data to use. Specifying TRUE will select all columns (default: All floating point value columns)
#' @param verbose        logical. If TRUE, show a progress bar and plot the output
#' 
#' @return Object of class \link{Sigmas}
#' 
#' @seealso \code{\link{Sigmas}}, the class returned by this; \code{\link{DiffusionMap}}, the class this is used for
#' 
#' @examples
#' data(guo)
#' sigs <- find.sigmas(guo, verbose = TRUE)
#' DiffusionMap(guo, sigs)
#' 
#' @export
find.sigmas <- function(
	data,
	step.size = 0.1,
	steps = 10L,
	start = NULL,
	sample.rows = 500L,
	early.exit = FALSE,
	...,
	censor.val = NULL, censor.range = NULL,
	missing.range = NULL,
	vars = NULL,
	verbose = TRUE
) {
	data <- extract.doublematrix(data, vars)
	
	if (any(is.na(data)))
		data <- as.matrix(hotdeck(data, imp_var = FALSE))
	
	stopifnot(steps >= 3L)
	
	if (length(sample.rows) > 1L) {
		data <- data[sample.rows, ]
	} else if (nrow(data) > sample.rows) {
		sample.idx <- sample(nrow(data), sample.rows)
		data <- data[sample.idx, ]
	}
	
	n <- nrow(data)
	
	dists <- dist(data)
	
	min.dist <- min(dists)
	if (min.dist == 0)
		stop('Minimum distance in the data may not be 0')
	
	if (is.null(start))
		start <- log10(min.dist)
	#if (missing(step.size))
	#  step.size = min.dist / steps
	
	if (verbose) print(c(min.dist = min.dist, start = start, step.size = step.size))
	
	get.trans.p <-
		if (test.censoring(censor.val, censor.range, data, missing.range)) {
			function(sigma) censoring(data, censor.val, censor.range, missing.range, sigma)
		} else {
			msqd <- as(-(dists ^ 2), 'Matrix')
			function(sigma) exp(msqd / (2 * sigma ^ 2))
		}
	
	do.step <- function(i) {
		# i can be negative!
		log.sigma <- start + i*step.size
		trans.p <- get.trans.p(10 ^ log.sigma)
		
		diag.d <- colSums(trans.p, na.rm = TRUE)
		
		list(avrd.norm = (sum(log10(diag.d/n) / diag.d)) / sum(1 / diag.d),
		     log.sigma = log.sigma)
	}
	
	avrd.norms <- numeric(steps)
	log.sigmas <- numeric(steps)
	dim.norms  <- numeric(steps - 1)
	
	step.diff <- function(step) {
		idxs <- c(step, step - 1)
		diff(avrd.norms[idxs]) / diff(log.sigmas[idxs])
	}
	
	a0 <- do.step(0L)
	avrd.norms[[1L]] <- a0$avrd.norm
	log.sigmas[[1L]] <- a0$log.sigma
	
	a1 <- do.step(1)
	dir <- 1L
	avrd.norms[[2L]] <- a1$avrd.norm
	log.sigmas[[2L]] <- a1$log.sigma
	
	if (step.diff(2L) < 0) {
		a1 <- do.step(-1L)
		dir <- -1L
		avrd.norms[[2L]] <- a1$avrd.norm
		log.sigmas[[2L]] <- a1$log.sigma
	}
	
	dim.norms[[1L]] <- step.diff(2L)
	
	for (step in seq(2L, steps)) {
		a.i = do.step(dir * (step - 1L))
		avrd.norms[[step]] <- a.i$avrd.norm
		log.sigmas[[step]] <- a.i$log.sigma
		
		dif.step <- step - 1
		dim.norms[[dif.step]] <- step.diff(step)
		
		if (early.exit && step > 2 && dim.norms[[dif.step]] < dim.norms[[dif.step - 1L]]) {
			avrd.norms <- avrd.norms[seq_len(step)]
			log.sigmas <- log.sigmas[seq_len(step)]
			dim.norms  <- dim.norms[seq_len(dif.step)]
			break
		}
	}
	
	if (early.exit && step == steps) warning('All steps were exhausted without finding a maximum. Using last encountered sigma')
	
	optimal.idx <- which.max(dim.norms)
	
	ret <- Sigmas(
		log.sigmas = log.sigmas,
		dim.norms  = dim.norms,
		optimal.sigma = 10 ^ mean(log.sigmas[c(optimal.idx, optimal.idx + 1L)]),
		optimal.idx   = optimal.idx,
		avrd.norms = avrd.norms)
	
	if (verbose) plot(ret)
	
	ret
}
