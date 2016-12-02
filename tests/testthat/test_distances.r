full_t_p_local <- function(dat, sigma, dists) {
	d2 <- as.matrix(dists ^ 2)
	
	S1 <- sigma %*% t(sigma)
	S2 <- outer(sigma ^ 2, sigma ^ 2, '+')
	
	rv <- sqrt(2 * S1 / S2) * exp(-d2 / (2*S2))
	rv[d2 == 0] <- 0
	rv
}

test_that('no_censoring produces the correct output for local sigma', {
	test_data <- matrix(rnorm(4L*5L), 4L, 5L)
	k <- 3L
	n_local <- 2:3
	
	knn <- find_knn(test_data, k)
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
