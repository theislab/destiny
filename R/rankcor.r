pr_rankcor <- function(x, y) {
	crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))
}

pr_rankcor_prefun <- function(x, y, pairwise, p, reg_entry) {
	x_ranked <- rank_mat(x)
	y_ranked <- rank_mat(y)
	x_rc <- sweep(x_ranked, 1, rowMeans(x_ranked))
	y_rc <- sweep(y_ranked, 1, rowMeans(y_ranked))
	list(x = x_rc, y = y_rc, pairwise = pairwise, p = p, reg_entry = reg_entry)
}

#' @importFrom proxy pr_DB
proxy_add_rankcor_simil <- function() {
	pr_DB$set_entry(
		FUN = pr_rankcor,
		PREFUN = pr_rankcor_prefun,
		names = c('rankcorrelation', 'rankcor', 'spearman'),
		distance = FALSE,
		convert = 'pr_simil2dist',
		type = 'metric',
		loop = TRUE,
		C_FUN = FALSE,
		abcd = FALSE,
		formula = 'xy / sqrt(xx * yy) for ranked and centered x,y',
		reference = '???',
		description = 'rank correlation (taking n instead of n-1 for the variance)')
}
