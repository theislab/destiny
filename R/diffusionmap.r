#' @import Matrix
#' @importFrom methods loadMethod
#' @importFrom BiocGenerics updateObject
#' @importFrom proxy dist
#' @importFrom FNN get.knn
#' @importFrom igraph arpack
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom VIM hotdeck
#' @useDynLib destiny
#' @include sigmas.r
#' @include s4-null-unions.r
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
#' @slot phi0           Density vector of normalized transition probability matrix
#' @slot d              Density vector of transition probability matrix
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
		phi0          = 'numeric',
		d             = 'numeric',
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
	
	sigmas <- sigma
	if (is.numeric(sigmas)) {
		sigmas <- new('Sigmas', 
			log.sigmas    = NULL,
			dim.norms     = NULL,
			optimal.sigma = sigma,
			optimal.idx   = NULL,
			avrd.norms    = NULL)
	} else if (!identical(distance, 'euclidean')) {
		stop(sprintf('You have to use euclidean distances with sigma estimation, not %s.', sQuote(distance)))
	} else if (is.null(sigmas)) {
		sigmas <- find.sigmas(
			imputed.data,
			distance = distance,
			censor.val = censor.val,
			censor.range = censor.range,
			missing.range = missing.range,
			vars = vars,
			verbose = verbose)
	} else if (!is(sigmas, 'Sigmas')) {
		stop(sprintf('The sigma parameter needs to be NULL, numeric or a %s object, not a %s.', sQuote('Sigmas'), sQuote(class(sigmas))))
	}
	sigma <- optimal.sigma(sigmas)
	
	# find KNNs
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
		trans.p <- censoring(data, censor.val, censor.range, missing.range, sigma, knn$nn.index, cb)
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
	rm(knn)  # free memory
	
	#nnzero
	
	# normalize trans.p and only retain intra-cell transitions
	diag(trans.p) <- 0
	trans.p <- drop0(trans.p)
	trans.p <- symmpart(trans.p) # double generic columnsparse to ... symmetric ...: dgCMatrix -> dsCMatrix
	
	d <- rowSums(trans.p, na.rm = TRUE)
	
	#ijk <- summary(trans.p)
	
	max.dist <- max(trans.p@x, na.rm = TRUE)
	stopifsmall(max.dist)
	
	if (density.norm) {
	  trans.p <- as(trans.p, 'dgTMatrix') # use non-symmetric triples to operate on all values
	  H <- sparseMatrix(trans.p@i, trans.p@j, x = trans.p@x / (d[trans.p@i + 1] * d[trans.p@j + 1]), dims = dim(trans.p), index1 = FALSE)
	  #creates a dgCMatrix
	} else {
	  H <- trans.p
	}
	rm(trans.p)  # free memory
	
	# only used for returning. could be used for (slower) eigen decomposition
	d_ <- rowSums(H)
	Hp <- Diagonal(x = d_ ^ -1) %*% H
	
	# calculate the inverse of a diagonal matrix by inverting the diagonal
	D.rot <- Diagonal(x = d_ ^ -.5)
	M <- D.rot %*% H %*% D.rot
	rm(H)  # free memory
	
	if (verbose) {
		cat('performing eigen decomposition...')
		tic <- proc.time()
	}
	#eig.M <- eig.decomp(Hp, n, n.eigs, FALSE)
	eig.M <- eig.decomp(M, n, n.eigs, TRUE)
	if (verbose) {
		cat('...done. Time:\n')
		print(proc.time() - tic)
	}
	
	#eig.vec <- eig.M$vectors
	eig.vec <- as.matrix(t(t(eig.M$vectors) %*% D.rot))
	colnames(eig.vec) <- paste0('DC', seq(0, n.eigs))
	
	new(
		'DiffusionMap',
		eigenvalues   = eig.M$values[-1],
		eigenvectors  = eig.vec[, -1],
		sigmas        = sigmas,
		data.env      = data.env,
		eigenvec0     = eig.vec[, 1],
		transitions   = Hp,
		phi0          = d_ / sum(d_),
		d             = d,
		k             = k,
		density.norm  = density.norm,
		distance      = distance,
		censor.val    = censor.val,
		censor.range  = censor.range,
		missing.range = missing.range)
}

extract.doublematrix <- function(data, vars = NULL) {
	if (is.data.frame(data)) {
		data <- as.matrix(data[sapply(data, is.double)])
	} else if (inherits(data, 'ExpressionSet')) {
		data <- t(exprs(data))
	} else if (!is.matrix(data)) {
		stop('Data needs to be matrix, data.frame or ExpressionSet')
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

eig.decomp <- function(M, n, n.eigs, sym) {
	f <- function(x, A = NULL) as.matrix(A %*% x)
	wh <- if (sym) 'LA' else 'LM'
	#constraints: n >= ncv > nev
	ar <- arpack(f, extra = M, sym = sym, options = list(
		which = wh, n = n, ncv = min(n, 4*n.eigs), nev = n.eigs + 1))
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
