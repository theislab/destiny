#' @include sigmas.r
#' @include s4-unions.r
#' @useDynLib destiny
NULL

sigma_msg <- function(sigma) sprintf(
	"The sigma parameter needs to be NULL, 'local', 'global', numeric or a %s object, not a %s.",
	sQuote('Sigmas'), sQuote(class(sigma)))

#' Create a diffusion map of cells
#' 
#' The provided data can be a double \link[base]{matrix} of expression data or a \link[base]{data.frame} with all non-integer (double) columns
#' being treated as expression data features (and the others ignored), an \link[Biobase:class.ExpressionSet]{ExpressionSet}, or a \link[SingleCellExperiment]{SingleCellExperiment}.
#' 
#' @param data           Expression data to be analyzed and covariates. Provide \code{vars} to select specific columns other than the default: all double value columns.
#'                       If \code{distance} is a distance matrix, \code{data} has to be a \code{\link{data.frame}} with covariates only.
#' @param sigma          Diffusion scale parameter of the Gaussian kernel. One of \code{'local'}, \code{'global'}, a (\link[base]{numeric}) global sigma or a \link{Sigmas} object.
#'                       When choosing \code{'global'}, a global sigma will be calculated using \code{\link{find_sigmas}}. (Optional. default: \code{'local'})
#'                       A larger sigma might be necessary if the eigenvalues can not be found because of a singularity in the matrix
#' @param k              Number of nearest neighbors to consider (default: a guess betweeen 100 and \eqn{n - 1}. See \code{\link{find_dm_k}}).
#' @param n_eigs         Number of eigenvectors/values to return (default: 20)
#' @param density_norm   logical. If TRUE, use density normalisation
#' @param ...            Unused. All parameters to the right of the \code{...} have to be specified by name (e.g. \code{DiffusionMap(data, distance = 'cosine')})
#' @param distance       Distance measurement method applied to \code{data} or a distance matrix/\code{\link[stats]{dist}}. For the allowed values, see \code{\link{find_knn}}.
#'                       If this is a \code{\link[Matrix]{sparseMatrix}}, zeros are interpreted as "not a close neighbors",
#'                       which allows the use of kNN-sparsified matrices (see the return value of \code{\link{find_knn}}.
#' @param n_local        If \code{sigma == 'local'}, the \code{n_local}th nearest neighbor(s) determine(s) the local sigma
#' @param rotate         logical. If TRUE, rotate the eigenvalues to get a slimmer diffusion map
#' @param censor_val     Value regarded as uncertain. Either a single value or one for every dimension (Optional, default: censor_val)
#' @param censor_range   Uncertainity range for censoring (Optional, default: none). A length-2-vector of certainty range start and end. TODO: also allow \eqn{2\times G} matrix
#' @param missing_range  Whole data range for missing value model. Has to be specified if NAs are in the data
#' @param vars           Variables (columns) of the data to use. Specifying NULL will select all columns (default: All floating point value columns)
#' @param verbose        Show a progressbar and other progress information (default: do it if censoring is enabled)
#' @param suppress_dpt   Specify TRUE to skip calculation of necessary (but spacious) information for \code{\link{DPT}} in the returned object (default: FALSE)
#' 
#' @return A DiffusionMap object:
#' 
#' @slot eigenvalues    Eigenvalues ranking the eigenvectors
#' @slot eigenvectors   Eigenvectors mapping the datapoints to \code{n_eigs} dimensions
#' @slot sigmas         \link{Sigmas} object with either information about the \link{find_sigmas} heuristic run or just local or \link{optimal_sigma}.
#' @slot data_env       Environment referencing the data used to create the diffusion map
#' @slot eigenvec0      First (constant) eigenvector not included as diffusion component.
#' @slot transitions    Transition probabilities. Can be NULL
#' @slot d              Density vector of transition probability matrix
#' @slot d_norm         Density vector of normalized transition probability matrix
#' @slot k              The k parameter for kNN
#' @slot n_local        The \code{n_local}th nearest neighbor(s) is/are used to determine local kernel density
#' @slot density_norm   Was density normalization used?
#' @slot rotate         Were the eigenvectors rotated?
#' @slot distance       Distance measurement method used
#' @slot censor_val     Censoring value
#' @slot censor_range   Censoring range
#' @slot missing_range  Whole data range for missing value model
#' @slot vars           Vars parameter used to extract the part of the data used for diffusion map creation
#' 
#' @seealso \link{DiffusionMap-methods} to get and set the slots. \code{\link{find_sigmas}} to pre-calculate a fitting global \code{sigma} parameter
#' 
#' @examples
#' data(guo)
#' DiffusionMap(guo)
#' DiffusionMap(guo, 13, censor_val = 15, censor_range = c(15, 40), verbose = TRUE)
#' 
#' covars <- data.frame(covar1 = letters[1:100])
#' dists <- dist(matrix(rnorm(100*10), 100))
#' DiffusionMap(covars, distance = dists)
#' 
#' @importFrom methods setClass validObject
#' @rdname DiffusionMap-class
#' @export
setClass(
	'DiffusionMap',
	slots = c(
		eigenvalues   = 'numeric',
		eigenvectors  = 'matrix',
		sigmas        = 'Sigmas',
		data_env      = 'environment',
		eigenvec0     = 'numeric',
		transitions   = 'dMatrixOrNULL',
		d             = 'numeric',
		d_norm        = 'numeric',
		k             = 'numeric',
		n_local       = 'numeric',
		density_norm  = 'logical',
		rotate        = 'logical',
		distance      = 'character',
		censor_val    = 'numericOrNULL',
		censor_range  = 'numericOrNULL',
		missing_range = 'numericOrNULL',
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
		#TODO: validate data_env (data, accumulated_transitions) and eigenvec0
		else if (length(object@d) != nrow(object@eigenvectors))
			'd must be as long as each eigenvector'
		else if (length(object@k) != 1)
			'k must be a number'
		else if (length(object@density_norm) != 1)
			'density_norm must be TRUE or FALSE'
		else if (length(object@rotate) != 1)
			'rotate must be TRUE or FALSE'
		else if (!(object@distance %in% c('euclidean', 'cosine', 'rankcor', 'l2', 'custom')))
			'distance must be "euclidean", "cosine", "rankcor", or "l2"'
		else if (is.null(object@censor_val) != is.null(object@censor_range))
			'Both censor_val and censor_range either need to be NULL or not'
		else if (!is.null(object@censor_val) && length(object@censor_val) != 1)
			'censor_val must be a number'
		else if (!is.null(object@censor_range) && (length(object@censor_range) != 2 || diff(object@censor_range) <= 0))
			'censor_range must be a increasing range (two numbers, the left one being larger)'
		else if (!is.null(object@missing_range) && (length(object@missing_range) != 2 || diff(object@missing_range) <= 0))
			'missing_range must be a increasing range (two numbers, the left one being larger)'
		#TODO validate vars
		else TRUE
	})

#' @importFrom methods new as is
#' @importFrom Matrix Diagonal colSums rowSums t
#' @importFrom VIM hotdeck
#' @rdname DiffusionMap-class
#' @export
DiffusionMap <- function(
	data = stopifnot_distmatrix(distance),
	sigma = 'local',
	k = find_dm_k(dataset_n_observations(data, distance) - 1L),
	n_eigs = min(20L, dataset_n_observations(data, distance) - 2L),
	density_norm = TRUE,
	...,
	distance = c('euclidean', 'cosine', 'rankcor', 'l2'),
	n_local = seq(to = min(k, 7L), length.out = min(k, 3L)),
	rotate = FALSE,
	censor_val = NULL, censor_range = NULL,
	missing_range = NULL,
	vars = NULL,
	verbose = !is.null(censor_range),
	suppress_dpt = FALSE
) {
	# make sure those promises are resolved before we mess with `data`
	force(k)
	force(n_eigs)
	
	stopifparams(...)
	
	if (is.null(sigma) || !is(sigma, 'Sigmas') && isTRUE(is.na(sigma)))
		sigma <- 'local'
	if (!is(sigma, 'Sigmas') && !(length(sigma) == 1L && sigma %in% c('local', 'global')) && !is.numeric(sigma))
		stop(sigma_msg(sigma))
	
	if (identical(sigma, 'local') && any(n_local > k))
		stop('For local sigma, All entries of n_local (', paste(n_local, collapse = ','), ') have to be \u2264 k (', k, ')')
	
	# store away data and continue using imputed, unified version
	data_env <- new.env(parent = .GlobalEnv)
	
	#TODO: SVD
	
	if (is_distmatrix(distance)) {
		if (!(is.data.frame(data) || is.null(data))) stop('If you provide a matrix for `distance`, `data` has to be NULL or a covariate `data.frame` is of class', class(data))
		
		data_env$data <- if (is.null(data)) distance else data  # put covariates or distance
		if (!is.null(rownames(data_env$data))) rownames(distance) <- colnames(distance) <- rownames(data)
		dists <- as(distance, 'symmetricMatrix')
		distance <- 'custom'
		imputed_data <- NULL
		n <- nrow(dists)
	} else {
		dists <- NULL
		distance <- match.arg(distance)
		
		data_env$data <- data
		data <- dataset_extract_doublematrix(data, vars)
		imputed_data <- data
		if (any(is.na(imputed_data)))
			imputed_data <- as.matrix(hotdeck(data, imp_var = FALSE))
		
		n <- nrow(imputed_data)
	}
	
	# arg validation
	
	if (n <= n_eigs + 1L) stop('Eigen decomposition not possible if n \u2264 n_eigs+1 (And ', n,' \u2264 ', n_eigs + 1L, ')')

	if (is.null(k) || is.na(k)) k <- n - 1L
	#TODO: optimize case
	#dense <- k == n - 1L
	
	if (k >= n) stop(sprintf('k has to be < nrow(data) (And %s \u2265 nrow(data))', k))
	
	censor <- test_censoring(censor_val, censor_range, imputed_data, missing_range)
	
	if (censor && !identical(distance, 'euclidean')) stop('censoring model only valid with euclidean distance')
	
	knn <- get_knn(imputed_data, dists, k, distance, verbose)  # use dists if given, else compute from data
	
	sigmas <- get_sigmas(imputed_data, knn$dist, sigma, n_local, distance, censor_val, censor_range, missing_range, vars, verbose)
	sigma <- optimal_sigma(sigmas)  # single number = global, multiple = local
	
	trans_p <- transition_probabilities(imputed_data, sigma, knn$dist_mat, censor, censor_val, censor_range, missing_range, verbose)
	rm(knn)  # free memory
	
	d <- rowSums(trans_p, na.rm = TRUE) + 1 # diagonal set to 1
	
	# normalize by density if requested
	norm_p <- get_norm_p(trans_p, d, d, density_norm)
	rm(trans_p)  # free memory
	
	d_norm <- rowSums(norm_p)
	
	# calculate the inverse of a diagonal matrix by inverting the diagonal
	d_rot <- Diagonal(x = d_norm ^ -.5)
	transitions <- as(d_rot %*% norm_p %*% d_rot, 'symmetricMatrix')
	rm(norm_p)  # free memory
	
	eig_transitions <- decomp_transitions(transitions, n_eigs + 1L, verbose)
	
	eig_vec <- eig_transitions$vectors
	eig_val <- eig_transitions$values
	if (rotate) eig_vec <- as.matrix(t(t(eig_vec) %*% d_rot))
	rownames(eig_vec) <-
		rownames(transitions) <-
		colnames(transitions) <-
		rownames(if (!is.null(rownames(dists))) dists else data)
	colnames(eig_vec) <-
		names(eig_val) <-
		paste0('DC', seq(0, n_eigs))
	
	new(
		'DiffusionMap',
		eigenvalues   = eig_val[-1],
		eigenvectors  = eig_vec[, -1, drop = FALSE],
		sigmas        = sigmas,
		data_env      = data_env,
		eigenvec0     = eig_vec[, 1],
		transitions   = if (suppress_dpt) NULL else transitions,
		d             = d,
		d_norm        = d_norm,
		k             = k,
		n_local       = n_local,
		rotate        = rotate,
		density_norm  = density_norm,
		distance      = distance,
		censor_val    = censor_val,
		censor_range  = censor_range,
		missing_range = missing_range)
}


is_distmatrix <- function(distance) {
	if (is.character(distance))
		FALSE
	else if ((is.matrix(distance) && diff(dim(distance)) == 0L) || is(distance, 'dist') || is(distance, 'symmetricMatrix'))
		TRUE
	else
		stop(
			'`distance` needs to be a distance measure or a square matrix, but is a ',
			class(dist), ' of dim ', nrow(distance), ' by ', ncol(distance))
}


stopifnot_distmatrix <- function(distance) {
	if (is_distmatrix(distance)) NULL
	else stop('`data`` needs to be set if `distance` is a distance measure')
}


#' @importFrom methods new is
get_sigmas <- function(imputed_data, nn_dists, sigma, n_local, distance = 'euclidean', censor_val = NULL, censor_range = NULL, missing_range = NULL, vars = NULL, verbose = FALSE) {
	unspecified_local <- identical(sigma, 'local')
	if (unspecified_local || is.numeric(sigma)) {
		if (unspecified_local) {
			sig_mat <- nn_dists[, n_local, drop = FALSE]
			sigma <- rowSums(sig_mat) / length(n_local) / 2
		}
		new('Sigmas', 
			log_sigmas    = NULL,
			dim_norms     = NULL,
			optimal_sigma = sigma,
			optimal_idx   = NULL,
			avrd_norms    = NULL)
	} else if (identical(sigma, 'global')) {
		if (!identical(distance, 'euclidean'))
			stop(sprintf('You have to use euclidean distances with sigma estimation, not %s.', sQuote(distance)))
		
		find_sigmas(
			imputed_data,
			distance = distance,
			censor_val = censor_val,
			censor_range = censor_range,
			missing_range = missing_range,
			vars = vars,
			verbose = verbose)
	} else if (is(sigma, 'Sigmas')) {
		sigma
	} else {
		stop(sigma_msg(sigma))
	}
}


get_knn <- function(imputed_data, dists, k, distance = 'euclidean', verbose = FALSE) {
	stopifnot(is.null(imputed_data) != is.null(dists))
	
	if (!is.null(dists)) {
		nn_dist <- t(apply(dists, 1, function(row) sort(row)[2:k]))
		list(dist = nn_dist, dist_mat = dists)
	} else {
		verbose_timing(verbose, 'finding knns', find_knn(imputed_data, k, distance = distance))
	}
}


#' @importFrom Matrix sparseMatrix diag<- drop0 forceSymmetric skewpart symmpart
#' @importFrom utils txtProgressBar setTxtProgressBar
transition_probabilities <- function(imputed_data, sigma, dists, censor, censor_val, censor_range, missing_range, verbose) {
	n <- nrow(dists)
	
	# create markovian transition probability matrix (trans_p)
	
	cb <- if (verbose) {
		pb <- txtProgressBar(1, n, style = 3)
		function(i) setTxtProgressBar(pb, i)
	} else invisible
	
	# initialize trans_p
	trans_p <- verbose_timing(verbose, 'Calculating transition probabilities', {
		if (censor)
			censoring(imputed_data, sigma, dists, censor_val, censor_range, missing_range, cb)
		else
			no_censoring(dists, sigma, cb)
	})
	
	if (verbose) close(pb)
	
	#nnzero
	
	# normalize trans_p and only retain intra-cell transitions
	diag(trans_p) <- 0
	trans_p <- drop0(trans_p)
	
	stopifnot(is(trans_p, 'symmetricMatrix'))
	trans_p
}

#' @importFrom Matrix sparseMatrix which tcrossprod
no_censoring <- function(dists, sigma, cb = invisible) {
	d2 <- dists ^ 2
	stopifnot(isSymmetric(d2))
	
	t_p <- if (length(sigma) == 1L) {
		exp(-d2@x / (2 * sigma ^ 2))
	} else {
		# TODO: optimize local sigma no-censoring case
		stopifnot(d2@uplo == 'U')
		mask <- d2 != 0 & upper.tri.sparse(dists)
		m <- function(mat) suppressMessages(as(mat, 'dsCMatrix')[mask])  # suppress warning about "inefficient .M.sub.i.logical"
		
		S1 <- m(tcrossprod(Matrix(sigma)))
		S2 <- m(outer(sigma ^ 2, sigma ^ 2, '+'))
		
		sqrt(2 * S1 / S2) * exp(-d2@x / S2)
	}
	
	sparseMatrix(d2@i, p = d2@p, x = t_p, dims = dim(d2), symmetric = TRUE, index1 = FALSE)
}


#' @importFrom methods as
#' @importFrom Matrix sparseMatrix
get_norm_p <- function(trans_p, d, d_new, density_norm) {
	if (density_norm) {
		trans_p <- as(trans_p, 'dgTMatrix') # use non-symmetric triples to operate on all values
		stopifsmall(max(trans_p@x, na.rm = TRUE))
		
		#creates a dgCMatrix
		sparseMatrix(trans_p@i, trans_p@j, x = trans_p@x / (d_new[trans_p@i + 1] * d[trans_p@j + 1]), dims = dim(trans_p), index1 = FALSE)
	} else {
		trans_p
	}
}


decomp_transitions <- function(transitions, n_eigs, verbose)
	verbose_timing(verbose, 'performing eigen decomposition', eig_decomp(transitions, n_eigs))
