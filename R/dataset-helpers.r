#' @importFrom Biobase exprs
#' @importFrom SummarizedExperiment assay colData
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
	} else if (inherits(data, 'SingleCellExperiment')) {
		#TODO: allow other name?
		data <- t(assay(data, 'logcounts'))
	} else if (!is.matrix(data)) {
		stop('Data needs to be matrix, data.frame, ExpressionSet, or SingleCellExperiment')
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


#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom pcaMethods pca scores
dataset_maybe_extract_pca <- function(data, n_pcs, verbose = FALSE) {
	stopifnot(is.null(n_pcs) || length(n_pcs) == 1L)
	# Suppress PCA computation
	if (!is.null(n_pcs) && is.na(n_pcs)) return(NULL)
	# If n_pcs is NULL, data needs to have a PCA
	has_pca <- inherits(data, 'SingleCellExperiment') && 'pca' %in% reducedDimNames(data)
	if (is.null(n_pcs) && !has_pca) return(NULL)
	# get PCs from SCE if possible
	if (has_pca) {
		pcs <- reducedDim(data, 'pca')
		if (is.null(n_pcs) || n_pcs == ncol(pcs)) {
			if (verbose) cat('Using reducedDim(data, "pca") to compute distances\n')
			return(pcs)
		} else if (n_pcs < ncol(pcs)) {
			warning('Specified n_pcs < ncol(reducedDim(data, "pca")), using subset')
			return(pcs[, seq_len(n_pcs), drop = FALSE])
		} else {# n_pcs > ncol(pcs)
			warning('Specified n_pcs > ncol(reducedDim(data, "pca")), recalculating PCA')
			rm(pcs)
		}
	}
	# No PCA in data or requested more PCs
	data_mat <- dataset_extract_doublematrix(data)
	if (ncol(data_mat) < n_pcs) stop('Cannot compute ', n_pcs, ' PCs from data with ', ncol(data_mat), ' columns.')
	scores(pca(data_mat, nPcs = n_pcs))
}


#' @importFrom Biobase sampleNames
dataset_n_observations <- function(data, distances) {
	if (is.null(data)) nrow(distances)
	else if (is(data, 'ExpressionSet')) length(sampleNames(data))
	else if (is(data, 'SingleCellExperiment')) ncol(data)
	else nrow(data)
}


#' @importFrom Biobase featureNames
dataset_n_features <- function(data, distances = NULL, vars = NULL) {
	if (is.null(data)) ncol(distances)
	else if (is(data, 'ExpressionSet')) length(featureNames(data))
	else if (is(data, 'SingleCellExperiment')) nrow(data)
	else if (!is.null(vars)) ncol(data[, vars])
	else if (is.data.frame(data)) ncol(data[, sapply(data, is.double)])
	else ncol(data)
}


#' @importFrom methods canCoerce
#' @importFrom utils getS3method
dataset_to_df <- function(dta, row.names = NULL, optional = FALSE, ...) {
	# The ExpressionSet as.data.frame sucks
	if (is(dta, 'ExpressionSet')) {
		cbind(as.data.frame(t(exprs(dta)), row.names, optional, ...), pData(dta))
	} else if (is(dta, 'SingleCellExperiment')) {
		smp_meta <- as.data.frame(colData(dta), row.names, optional, ...)
		
		#TODO: allow other name?
		mat <- assay(dta, 'logcounts')
		if (is(mat, 'sparseMatrix')) {
			n_genes_max <- getOption('destiny.genes.max', 5000L)
			if (nrow(mat) > n_genes_max) {
				warning(
					'Data has too many genes to convert it to a data.frame (',
					nrow(mat), ' > ', n_genes_max,
					'). You can try setting options(destiny.genes.max = ...)')
				return(smp_meta)
			} else {
				mat <- as.matrix(mat)
			}
		}
		
		cbind(as.data.frame(t(mat), row.names, optional, ...), smp_meta)
	} else if (canCoerce(dta, 'data.frame')) {
		as(dta, 'data.frame')
	} else if (!is.null(getS3method('as.data.frame', class(dta)[[1L]], optional = TRUE))) {
		as.data.frame(dta, row.names, optional, ...)
	} else NULL
}


dataset_names <- function(dta) {
	if (is(dta, 'ExpressionSet')) {
		c(featureNames(dta), varLabels(dta))
	} else if (is(dta, 'SingleCellExperiment')) {
		c(rownames(dta), colnames(colData(dta)))
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
		} else if (is(dta, 'SingleCellExperiment') && f %in% rownames(dta)) {
			#TODO: allow other name?
			assay(dta, 'logcounts')[f, ]
		} else if (is(dta, 'ExpressionSet') || is(dta, 'SingleCellExperiment') || is.data.frame(dta)) {
			# data.frame or phenoData/colData
			dta[[f]]
		} else {
			dta[, f]
		}
	}, error = function(e) stop(sprintf('Invalid `%s`: No column, annotation, or feature found with name %s', deparse(substitute(f)), dQuote(f))))
}
