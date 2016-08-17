#' @include sigmas.r
#' @include s4-null-unions.r
#' @useDynLib destiny
NULL

sigma_msg <- function(sigma) sprintf(
	"The sigma parameter needs to be NULL, 'local', 'global', numeric or a %s object, not a %s.",
	sQuote('Sigmas'), sQuote(class(sigma)))

#' Create a diffusion map of cells
#' 
#' The provided data can be a double \link[base]{matrix} of expression data or a \link[base]{data.frame} with all non-integer (double) columns
#' being treated as expression data features (and the others ignored), or an \link[Biobase]{ExpressionSet}.
#' 
#' @param data           Expression data to be analyzed. Provide \code{vars} to select specific columns other than the default: all double value columns
#' @param sigma          Diffusion scale parameter of the Gaussian kernel. One of \code{'local'}, \code{'global'}, a (\link[base]{numeric}) global sigma or a \link{Sigmas} object.
#'                       When choosing \code{'global'}, a global sigma will be calculated using \code{\link{find_sigmas}}. (Optional. default: \code{'local'})
#'                       A larger sigma might be necessary if the eigenvalues can not be found because of a singularity in the matrix
#' @param k              Number of nearest neighbors to consider (default: a guess betweeen 100 and \eqn{n - 1}. See \code{\link{find_dm_k}}).
#' @param n_eigs         Number of eigenvectors/values to return (default: 20)
#' @param density_norm   logical. If TRUE, use density normalisation
#' @param ...            All parameter after this are optional and have to be specified by name
#' @param distance       Distance measurement method. Euclidean distance (default), cosine distance (\eqn{1-corr(c_1, c_2)}) or rank correlation distance (\eqn{1-corr(rank(c_1), rank(c_2))}).
#' @param n_local        If \code{sigma == 'local'}, the \code{n_local}th nearest neighbor determines the local sigma.
#' @param censor_val     Value regarded as uncertain. Either a single value or one for every dimension (Optional, default: censor_val)
#' @param censor_range   Uncertainity range for censoring (Optional, default: none). A length-2-vector of certainty range start and end. TODO: also allow \eqn{2\times G} matrix
#' @param missing_range  Whole data range for missing value model. Has to be specified if NAs are in the data
#' @param vars           Variables (columns) of the data to use. Specifying NULL will select all columns (default: All floating point value columns)
#' @param verbose        Show a progressbar and other progress information (default: do it if censoring is enabled)
#' @param suppress_dpt   Specify TRUE to skip calculation of necessary (but spacious) information for \code{\link{dpt}} in the returned object (default: FALSE)
#' 
#' @return A DiffusionMap object:
#' 
#' @slot eigenvalues    Eigenvalues ranking the eigenvectors
#' @slot eigenvectors   Eigenvectors mapping the datapoints to \code{n_eigs} dimensions
#' @slot sigmas         \link{Sigmas} object with either information about the \link{find_sigmas} heuristic run or just local or \link{optimal_sigma}.
#' @slot data_env       Environment referencing the data used to create the diffusion map
#' @slot eigenvec0      First (constant) eigenvector not included as diffusion component.
#' @slot transitions    Transition probabilities. Can be NULL
#' @slot propagations   Propagation matrix derived from transition probabilites. Can be NULL
#' @slot d              Density vector of transition probability matrix
#' @slot d_norm         Density vector of normalized transition probability matrix
#' @slot k              The k parameter for kNN
#' @slot density_norm   Was density normalization used?
#' @slot distance       Distance measurement method used.
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
		data_env      = 'environment',
		eigenvec0     = 'numeric',
		transitions   = 'dMatrixOrNULL',
		propagations  = 'dMatrixOrNULL',
		d             = 'numeric',
		d_norm        = 'numeric',
		k             = 'numeric',
		n_local       = 'numeric',
		density_norm  = 'logical',
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
		#TODO: validate data_env and eigenvec0
		else if (length(object@d) != nrow(object@eigenvectors))
			'd must be as long as each eigenvector'
		else if (length(object@k) != 1)
			'k must be a number'
		else if (length(object@density_norm) != 1)
			'density_norm must be TRUE or FALSE'
		else if (!(object@distance %in% c('euclidean', 'cosine', 'rankcor')))
			'distance must be "euclidean", "cosine" or "rankcor"'
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
#' @name DiffusionMap class
#' @export
DiffusionMap <- function(
	data,
	sigma = 'local',
	k = find_dm_k(nrow(data) - 1L),
	n_eigs = min(20L, nrow(data) - 2L),
	density_norm = TRUE,
	...,
	distance = c('euclidean', 'cosine', 'rankcor'),
	n_local = 5L,
	censor_val = NULL, censor_range = NULL,
	missing_range = NULL,
	vars = NULL,
	verbose = !is.null(censor_range),
	suppress_dpt = FALSE
) {
	if (is.null(sigma) || isTRUE(is.na(sigma)))
		sigma <- 'local'
	if (!(length(sigma) == 1L && sigma %in% c('local', 'global')) && !is.numeric(sigma) && !is(sigma, 'Sigmas'))
		stop(sigma_msg(sigma))
	
	distance <- match.arg(distance)
	
	# store away data and continue using imputed, unified version
	data_env <- new.env(parent = .GlobalEnv)
	data_env$data <- data
	
	data <- extract_doublematrix(data, vars)
	
	#TODO: SVD
	
	imputed_data <- data
	if (any(is.na(imputed_data)))
		imputed_data <- as.matrix(hotdeck(data, imp_var = FALSE))
	
	n <- nrow(imputed_data)
	
	# arg validation
	
	if (n <= n_eigs + 1L) stop(sprintf('Eigen decomposition not possible if n \u2264 n_eigs+1 (And %s \u2264 %s)', n, n_eigs + 1L))
	
	if (is.null(k) || is.na(k)) k <- n - 1L
	#TODO: optimize case
	#dense <- k == n - 1L
	
	if (k >= nrow(imputed_data)) stop(sprintf('k has to be < nrow(data) (And %s \u2265 nrow(data))', k))
	
	censor <- test_censoring(censor_val, censor_range, data, missing_range)
	
	if (censor && !identical(distance, 'euclidean')) stop('censoring model only valid with euclidean distance')
	
	knn <- find_knn(imputed_data, k, verbose)
	
	sigmas <- get_sigmas(imputed_data, sigma, distance, knn, n_local, censor_val, censor_range, missing_range, vars, verbose)
	sigma <- optimal_sigma(sigmas)  # single number = global, multiple = local
	
	trans_p <- transition_probabilities(imputed_data, sigma, distance, knn, censor, censor_val, censor_range, missing_range, verbose)
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
	
	eig_transitions <- decomp_transitions(transitions, n_eigs, verbose)
	
	eig_vec <- as.matrix(t(t(eig_transitions$vectors) %*% d_rot))
	colnames(eig_vec) <- paste0('DC', seq(0, n_eigs))
	
	if (suppress_dpt) transitions <- NULL
	propagations <- get_propagation_matrix(transitions, d_norm)
	
	new(
		'DiffusionMap',
		eigenvalues   = eig_transitions$values[-1],
		eigenvectors  = eig_vec[, -1],
		sigmas        = sigmas,
		data_env      = data_env,
		eigenvec0     = eig_vec[, 1],
		transitions   = transitions,
		propagations  = propagations,
		d             = d,
		d_norm        = d_norm,
		k             = k,
		n_local       = n_local,
		density_norm  = density_norm,
		distance      = distance,
		censor_val    = censor_val,
		censor_range  = censor_range,
		missing_range = missing_range)
}


#' @importFrom methods new is
get_sigmas <- function(imputed_data, sigma, distance, knn, n_local, censor_val, censor_range, missing_range, vars, verbose) {
	unspecified_local <- identical(sigma, 'local')
	if (unspecified_local || is.numeric(sigma)) {
		if (unspecified_local)
			sigma <- knn$nn.dist[, n_local] / 2
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


# get.knn(...)$nn.index : \eqn{n \times k} matrix for the nearest neighbor indices
# get.knn(...)$nn.dist  : \eqn{n \times k} matrix for the nearest neighbor Euclidean distances.
#' @importFrom FNN get.knn
find_knn <- function(imputed_data, k, verbose)
	verbose_timing(verbose, 'finding knns', get.knn(imputed_data, k, algorithm = 'cover_tree'))


#' @importFrom Matrix sparseMatrix diag<- drop0 forceSymmetric skewpart symmpart
#' @importFrom utils txtProgressBar setTxtProgressBar
transition_probabilities <- function(imputed_data, sigma, distance, knn, censor, censor_val, censor_range, missing_range, verbose) {
	n <- nrow(knn$nn.index)
	
	# create markovian transition probability matrix (trans_p)
	
	cb <- if (verbose) {
		pb <- txtProgressBar(1, n, style = 3)
		function(i) setTxtProgressBar(pb, i)
	} else invisible
	
	# initialize trans_p
	trans_p <- verbose_timing(verbose, 'Calculating transition probabilities', {
		if (censor)
			censoring(imputed_data, censor_val, censor_range, missing_range, sigma, knn$nn.index, cb)
		else
			no_censoring(imputed_data, sigma, distance, knn, cb)
	})
	
	if (verbose) close(pb)
	
	#nnzero
	
	# normalize trans_p and only retain intra-cell transitions
	diag(trans_p) <- 0
	trans_p <- drop0(trans_p)
	
	# (double generic columnsparse to ... symmetric ...: dgCMatrix -> dsCMatrix)
	# retain all differences fully. symmpart halves them in the case of trans_p[i,j] == 0 && trans_p[j,i] > 0
	trans_p <- symmpart(trans_p) + abs(forceSymmetric(skewpart(trans_p))) # TODO: more efficient
	
	trans_p
}

#' @importFrom Matrix sparseMatrix which
no_censoring <- function(imputed_data, sigma, distance, knn, cb) {
	d2 <- switch(distance,
		euclidean  = d2_no_censor(knn$nn.index, knn$nn.dist,  cb),
		cosine  = icor2_no_censor(knn$nn.index, imputed_data, cb),
		rankcor = icor2_no_censor(knn$nn.index, imputed_data, cb, TRUE))
	
	t_p <- if (length(sigma) == 1L) {
		exp(-d2@x / (2 * sigma ^ 2))
	} else {
		# TODO: optimize local sigma no-censoring case
		mask <- which(d2 == 0)
		S1 <- (sigma %*% t(sigma))[mask]
		S2 <- outer(sigma ^ 2, sigma ^ 2, '+')[mask]
		sqrt(2 * S1 / S2) * exp(-d2@x / S2)
	}
	
	sparseMatrix(d2@i, p = d2@p, x = t_p, dims = dim(d2), index1 = FALSE)
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


#' @importFrom Matrix Diagonal
#' @importMethodsFrom Matrix solve
get_propagation_matrix <- function(transitions, d_norm) {
	if (is.null(transitions)) return(NULL)
	
	phi0 <- d_norm / sqrt(sum(d_norm ^ 2))
	d_rot <- Diagonal(x = d_norm ^ -.5)
	inv <- solve(Diagonal(n) - transitions + phi0 %*% t(phi0))
	inv - Diagonal(n)
}
