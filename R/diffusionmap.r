#' @include sigmas.r
#' @include s4-null-unions.r
#' @useDynLib destiny
NULL

#' Create a diffusion map of cells
#' 
#' The provided data can be a double \link[base]{matrix} of expression data or a \link[base]{data.frame} with all non-integer (double) columns
#' being treated as expression data features (and the others ignored), or an \link[Biobase]{ExpressionSet}.
#' 
#' @param data           Expression data to be analyzed. Provide \code{vars} to select specific columns other than the default: all double value columns
#' @param sigma          Diffusion scale parameter of the Gaussian kernel. Either a number or a \link{Sigmas} object.
#'                       (Optional. default: will be calculated using \link{find.sigmas})
#'                       A larger sigma might be necessary if the eigenvalues can not be found because of a singularity in the matrix
#' @param k              Number of nearest neighbors to consider (default: a guess betweeen 100 and \eqn{n - 1})
#'                       \code{NULL} or \code{NA} are also interpreted as \code{n - 1L}.
#' @param n.eigs         Number of eigenvectors/values to return (default: 20)
#' @param density.norm   logical. If TRUE, use density normalisation
#' @param ...            All parameter after this are optional and have to be specified by name
#' @param distance       Distance measurement method. Euclidean distance (default), cosine distance (\eqn{1-corr(c_1, c_2)}) or rank correlation distance (\eqn{1-corr(rank(c_1), rank(c_2))}).
#' @param censor.val     Value regarded as uncertain. Either a single value or one for every dimension (Optional, default: CENSOR.VAL)
#' @param censor.range   Uncertainity range for censoring (Optional, default: none). A length-2-vector of certainty range start and end. TODO: also allow \eqn{2\times G} matrix
#' @param missing.range  Whole data range for missing value model. Has to be specified if NAs are in the data
#' @param vars           Variables (columns) of the data to use. Specifying NULL will select all columns (default: All floating point value columns)
#' @param verbose        Show a progressbar and other progress information (default: do it if censoring is enabled)
#' 
#' @return A DiffusionMap object:
#' 
#' @slot eigenvalues    Eigenvalues ranking the eigenvectors
#' @slot eigenvectors   Eigenvectors mapping the datapoints to \code{n.eigs} dimensions
#' @slot sigmas         \link{Sigmas} object with either information about the \link{find.sigmas} heuristic run or just \link{optimal.sigma}.
#' @slot data.env       Environment referencing the data used to create the diffusion map
#' @slot eigenvec0      First (constant) eigenvector not included as diffusion component.
#' @slot transitions    Transition probabilities
#' @slot d              Density vector of transition probability matrix
#' @slot d_norm         Density vector of normalized transition probability matrix
#' @slot k              The k parameter for kNN
#' @slot density.norm   Was density normalization used?
#' @slot distance       Distance measurement method used.
#' @slot censor.val     Censoring value
#' @slot censor.range   Censoring range
#' @slot missing.range  Whole data range for missing value model
#' @slot vars           Vars parameter used to extract the part of the data used for diffusion map creation
#' 
#' @seealso \link{DiffusionMap-methods} to get and set the slots. \code{\link{find.sigmas}} to pre-calculate a fitting \code{sigma} parameter
#' 
#' @examples
#' data(guo)
#' DiffusionMap(guo)
#' DiffusionMap(guo, 13, censor.val = 15, censor.range = c(15, 40), verbose = TRUE)
#' 
#' @aliases DiffusionMap DiffusionMap-class
#' 
#' @importFrom methods setClass validObject
#' @name DiffusionMap class
#' @export
setClass(
	'DiffusionMap',
	slots = c(
		eigenvalues   = 'numeric',
		eigenvectors  = 'matrix',
		sigmas        = 'Sigmas',
		data.env      = 'environment',
		eigenvec0     = 'numeric',
		transitions   = 'dMatrix',
		d             = 'numeric',
		d_norm        = 'numeric',
		k             = 'numeric',
		density.norm  = 'logical',
		distance      = 'character',
		censor.val    = 'numericOrNULL',
		censor.range  = 'numericOrNULL',
		missing.range = 'numericOrNULL',
		vars          = 'characterOrnumericOrNULL'),
	validity = function(object) {
		# don't use getters as thisfunction is used in them!
		if (!is.vector(object@eigenvalues) || !is.numeric(object@eigenvalues))
			'eigenvalues(dm) has to be a numeric vector'
		else if (!is.matrix(object@eigenvectors) || !is.numeric(object@eigenvectors))
			'eigenvectors(dm) has to be a numeric matrix'
		else if (length(object@eigenvalues) != ncol(object@eigenvectors))
			'There must be exactly one eigenvalue per eigenvector: length(eigenvalues(dm)) == ncol(eigenvectors(dm))'
		else if (!isTRUE(validObject(object@sigmas, test = TRUE)))
			paste('sigmas invalid:', validObject(object@sigmas, test = TRUE))
		#TODO: validate data.env and eigenvec0
		else if (length(object@d) != nrow(object@eigenvectors))
			'd must be as long as each eigenvector'
		else if (length(object@k) != 1)
			'k must be a number'
		else if (length(object@density.norm) != 1)
			'density.norm must be TRUE or FALSE'
		else if (!(object@distance %in% c('euclidean', 'cosine', 'rankcor')))
			'distance must be "euclidean", "cosine" or "rankcor"'
		else if (is.null(object@censor.val) != is.null(object@censor.range))
			'Both censor.val and censor.range either need to be NULL or not'
		else if (!is.null(object@censor.val) && length(object@censor.val) != 1)
			'censor.val must be a number'
		else if (!is.null(object@censor.range) && (length(object@censor.range) != 2 || diff(object@censor.range) <= 0))
			'censor.range must be a increasing range (two numbers, the left one being larger)'
		else if (!is.null(object@missing.range) && (length(object@missing.range) != 2 || diff(object@missing.range) <= 0))
			'missing.range must be a increasing range (two numbers, the left one being larger)'
		#TODO validate vars
		else TRUE
	})

#' @importFrom methods new as
#' @importFrom Matrix Diagonal colSums rowSums t
#' @importFrom VIM hotdeck
#' 
#' @name DiffusionMap class
#' @export
DiffusionMap <- function(
	data,
	sigma = NULL,
	k = find.dm.k(nrow(data) - 1L),
	n.eigs = min(20L, nrow(data) - 2L),
	density.norm = TRUE,
	...,
	distance = c('euclidean', 'cosine', 'rankcor'),
	censor.val = NULL, censor.range = NULL,
	missing.range = NULL,
	vars = NULL,
	verbose = !is.null(censor.range)
) {
	distance <- match.arg(distance)
	
	# store away data and continue using imputed, unified version
	data.env <- new.env(parent = .GlobalEnv)
	data.env$data <- data
	
	data <- extract.doublematrix(data, vars)
	
	#TODO: SVD
	
	imputed.data <- data
	if (any(is.na(imputed.data)))
		imputed.data <- as.matrix(hotdeck(data, imp_var = FALSE))
	
	n <- nrow(imputed.data)
	
	# arg validation
	
	if (n <= n.eigs + 1L) stop(sprintf('Eigen decomposition not possible if n \u2264 n.eigs+1 (And %s \u2264 %s)', n, n.eigs + 1L))
	
	if (is.null(k) || is.na(k)) k <- n - 1L
	#TODO: optimize case
	#dense <- k == n - 1L
	
	if (k >= nrow(imputed.data)) stop(sprintf('k has to be < nrow(data) (And %s \u2265 nrow(data))', k))
	
	censor <- test.censoring(censor.val, censor.range, data, missing.range)
	
	if (censor && !identical(distance, 'euclidean')) stop('censoring model only valid with euclidean distance') 
	
	sigmas <- get.sigmas(imputed.data, sigma, distance, censor.val, censor.range, missing.range, vars, verbose)
	sigma <- optimal.sigma(sigmas)
	
	knn <- find.knn(imputed.data, k, verbose)
	
	trans.p <- transition.probabilities(imputed.data, distance, sigma, knn, censor, censor.val, censor.range, missing.range, verbose)
	rm(knn)  # free memory
	
	d <- rowSums(trans.p, na.rm = TRUE) + 1 # diagonal set to 1
	
	# normalize by density if requested
	norm_p <- get_norm_p(trans.p, d, d, density.norm)
	rm(trans.p)  # free memory
	
	d_norm <- rowSums(norm_p)
	
	# calculate the inverse of a diagonal matrix by inverting the diagonal
	D.rot <- Diagonal(x = d_norm ^ -.5)
	transitions <- as(D.rot %*% norm_p %*% D.rot, 'symmetricMatrix')
	rm(norm_p)  # free memory
	
	eig_transitions <- decomp_transitions(transitions, n.eigs, verbose)
	
	eig_vec <- as.matrix(t(t(eig_transitions$vectors) %*% D.rot))
	colnames(eig_vec) <- paste0('DC', seq(0, n.eigs))
	
	new(
		'DiffusionMap',
		eigenvalues   = eig_transitions$values[-1],
		eigenvectors  = eig_vec[, -1],
		sigmas        = sigmas,
		data.env      = data.env,
		eigenvec0     = eig_vec[, 1],
		transitions   = transitions,
		d             = d,
		d_norm        = d_norm,
		k             = k,
		density.norm  = density.norm,
		distance      = distance,
		censor.val    = censor.val,
		censor.range  = censor.range,
		missing.range = missing.range)
}

#' @importFrom Biobase exprs
extract.doublematrix <- function(data, vars = NULL) {
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

stopifsmall <- function(max.dist) {
	if (max.dist < .Machine$double.eps)
		stop(sprintf(
			'The supplied sigma is not large enough. Please select a larger one.
			find.sigmas(data) should return one with the right order of magnitude. (max dist. is %.3e)',
			max.dist))
}


#' Fast eigen decomposition using ARPACK
#'
#' @param M       A matrix (e.g. from the Matrix package)
#' @param n_eigs  Number of eigenvectors to return
#' @param sym     TRUE if M is symmetric
#' 
#' @return n eigenvectors of the transition matrix
#' 
#' @importFrom Matrix isSymmetric
#' @importFrom igraph arpack
#' @export
eig_decomp <- function(M, n_eigs, sym = isSymmetric(M)) {
	n <- nrow(M)
	f <- function(x, A = NULL) as.matrix(A %*% x)
	wh <- if (sym) 'LA' else 'LM'
	#constraints: n >= ncv > nev
	ar <- arpack(f, extra = M, sym = sym, options = list(
		which = wh, n = n, ncv = min(n, 4*n_eigs), nev = n_eigs + 1))
	if (!sym) {
		ar$vectors <- Re(ar$vectors)
		ar$values  <- Re(ar$values)
	}
	ar
}

#' Find a suitable k
#' 
#' The \code{k} parameter for the k nearest neighbors used in \link{DiffusionMap} should be as big as possible while
#' still being computationally feasible. This function approximates it depending on the size of the dataset \code{n}.
#' 
#' @param n      Number of possible neighbors (nrow(dataset) - 1)
#' @param min.k  Minimum number of neighbors. Will be chosen for \eqn{n \ge big}
#' @param small  Number of neighbors considered small. If/where \eqn{n \le small}, n itself will be returned.
#' @param big    Number of neighbors considered big. If/where \eqn{n \ge big}, \code{min.k} will be returned.
#' 
#' @return A vector of the same length as \code{n} that contains suitable \code{k} values for the respective \code{n}
#' 
#' @examples
#' curve(find.dm.k(n),     0, 13000, xname = 'n')
#' curve(find.dm.k(n) / n, 0, 13000, xname = 'n')
#' @export
find.dm.k <- function(n, min.k = 100L, small = 1000L, big = 10000L) {
	stopifnot(small < big)
	if (is.null(n)) return(NULL)
	
	k <- rep(NA_integer_, length(n))
	k[small >= n] <- n[small >= n]
	k[n >= big]   <- min.k
	
	rest <- !is.na(n) & small < n & n < big
	
	n.shifted <- (n[rest] - small) / (big - small)     # linear transf [small, big] -> [0, 1]
	k.shifted <- (cos(n.shifted * pi) + 1) / 2         # ease function [0, 1] -> [1, 0]
	k.rest    <- min.k + k.shifted * (n[rest] - min.k) # linear transf [0, 1] -> [min.k, n]
	
	k[rest] <- as.integer(round(k.rest))
	
	k
}

#' @importFrom methods new is
get.sigmas <- function(imputed.data, sigma, distance, censor.val, censor.range, missing.range, vars, verbose) {
	sigmas <- sigma
	if (is.numeric(sigmas)) {
		new('Sigmas', 
			log.sigmas    = NULL,
			dim.norms     = NULL,
			optimal.sigma = sigma,
			optimal.idx   = NULL,
			avrd.norms    = NULL)
	} else if (!identical(distance, 'euclidean')) {
		stop(sprintf('You have to use euclidean distances with sigma estimation, not %s.', sQuote(distance)))
	} else if (is.null(sigmas)) {
		find.sigmas(
			imputed.data,
			distance = distance,
			censor.val = censor.val,
			censor.range = censor.range,
			missing.range = missing.range,
			vars = vars,
			verbose = verbose)
	} else if (is(sigmas, 'Sigmas')) {
		sigmas
	} else {
		stop(sprintf('The sigma parameter needs to be NULL, numeric or a %s object, not a %s.', sQuote('Sigmas'), sQuote(class(sigmas))))
	}
}

#' @importFrom FNN get.knn
find.knn <- function(imputed.data, k, verbose) {
	if (verbose) {
		cat('finding knns...')
		tic <- proc.time()
	}
	knn <- get.knn(imputed.data, k, algorithm = 'cover_tree')
	# knn$nn.index : \eqn{n \times k} matrix for the nearest neighbor indices
	# knn$nn.dist  : \eqn{n \times k} matrix for the nearest neighbor Euclidean distances.
	if (verbose) {
		cat('...done. Time:\n')
		print(proc.time() - tic)
	}
	
	knn
}

#' @importFrom Matrix sparseMatrix drop0 forceSymmetric skewpart symmpart
#' @importFrom utils txtProgressBar setTxtProgressBar
transition.probabilities <- function(imputed.data, distance, sigma, knn, censor, censor.val, censor.range, missing.range, verbose) {
	n <- nrow(knn$nn.index)
	
	# create markovian transition probability matrix (trans.p)
	cb <- invisible
	if (verbose) {
		pb <- txtProgressBar(1, n, style = 3)
		cb <- function(i) setTxtProgressBar(pb, i)
		cat('Calculating transition probabilities...\n')
		tic <- proc.time()
	}
	
	# initialize trans.p
	if (censor) {
		trans.p <- censoring(imputed.data, censor.val, censor.range, missing.range, sigma, knn$nn.index, cb)
	} else {
		d2 <- switch(distance,
			euclidean  = d2_no_censor(knn$nn.index, knn$nn.dist,  cb),
			cosine  = icor2_no_censor(knn$nn.index, imputed.data, cb),
			rankcor = icor2_no_censor(knn$nn.index, imputed.data, cb, TRUE))
		
		trans.p <- sparseMatrix(d2@i, p = d2@p, x = exp(-d2@x / (2 * sigma ^ 2)), dims = dim(d2), index1 = FALSE)
		rm(d2)
	}
	
	if (verbose) {
		close(pb)
		cat('...done. Time:\n')
		print(proc.time() - tic)
	}
	
	#nnzero
	
	# normalize trans.p and only retain intra-cell transitions
	diag(trans.p) <- 0
	trans.p <- drop0(trans.p)
	
	# (double generic columnsparse to ... symmetric ...: dgCMatrix -> dsCMatrix)
	# retain all differences fully. symmpart halves them in the case of trans.p[i,j] == 0 && trans.p[j,i] > 0
	trans.p <- symmpart(trans.p) + abs(forceSymmetric(skewpart(trans.p))) # TODO: more efficient
	
	trans.p
}

#' @importFrom methods as
#' @importFrom Matrix sparseMatrix
get_norm_p <- function(trans.p, d, d_new, density.norm) {
	if (density.norm) {
		trans.p <- as(trans.p, 'dgTMatrix') # use non-symmetric triples to operate on all values
		stopifsmall(max(trans.p@x, na.rm = TRUE))
		
		#creates a dgCMatrix
		sparseMatrix(trans.p@i, trans.p@j, x = trans.p@x / (d_new[trans.p@i + 1] * d[trans.p@j + 1]), dims = dim(trans.p), index1 = FALSE)
	} else {
		trans.p
	}
}

decomp_transitions <- function(transitions, n_eigs, verbose) {
	if (verbose) {
		cat('performing eigen decomposition...')
		tic <- proc.time()
	}
	eig_transitions <- eig_decomp(transitions, n_eigs)
	if (verbose) {
		cat('...done. Time:\n')
		print(proc.time() - tic)
	}
	
	eig_transitions
}
