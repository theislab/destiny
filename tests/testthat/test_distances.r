library(Matrix)

full_t_p_local <- function(dat, sigma, dists) {
	d2 <- as.matrix(dists ^ 2)
	
	S1 <- sigma %*% t(sigma)
	S2 <- outer(sigma ^ 2, sigma ^ 2, '+')
	
	rv <- sqrt(2 * S1 / S2) * exp(-d2 / S2)
	rv[d2 == 0] <- 0
	rv
}

test_that('no_censoring produces the correct output for local sigma', {
	test_data <- matrix(rnorm(4L*5L), 4L, 5L)
	k <- 3L
	n_local <- 2:3
	
	knn <- find_knn(test_data, NULL, k)
	expect_identical(dim(knn$nn_dist), c(nrow(test_data), k))
	expect_identical(dim(knn$dist), rep(nrow(test_data), 2L))
	expect_identical(sum(!is.finite(knn$dist)), 0L)
	
	sigma <- optimal_sigma(get_sigmas(test_data, knn$nn_dist, 'local', n_local))
	expect_identical(length(sigma), nrow(test_data))
	
	dists_expected <- full_t_p_local(test_data, sigma, knn$dist)
	dists <- no_censoring(test_data, sigma, knn$dist)
	expect_equal(dim(dists), dim(dists_expected))
	expect_equal(as.matrix(dists), dists_expected)
})

matidx_apply <- function(nrow, ncol, FUN) {
	mat <- matrix(NA, nrow, ncol)
	vals <- mapply(FUN, row(mat), col(mat))
	matrix(vals, nrow, ncol)
}

test_that('no_censoring produces the correct output for cor distance', {
	test_data <- matrix(rnorm(4L*5L), 4L, 5L)
	k <- 3L
	
	cor_dists_expected <- matidx_apply(nrow(test_data), nrow(test_data), function(r1, r2) 1 - cor(test_data[r1, ], test_data[r2, ]))
	cor_dists_expected <- cor_dists_expected ^ 2
	
	dists <- destiny:::find_knn(test_data, NULL, k)$dist
	cor_dists <- destiny:::no_censoring_icor_rank(dists, test_data)
	cor_dists <- as.matrix(cor_dists)
	dimnames(cor_dists) <- NULL
	
	expect_identical(dim(cor_dists), dim(cor_dists_expected))
	expect_equal(cor_dists, cor_dists_expected)
})
