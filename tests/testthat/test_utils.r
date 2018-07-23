context('utils')

library(Matrix)

test_that('duplicated.dgCMatrix works', {
	m <- Matrix(c(
		1, 1, 2,
		0, 0, 0,
		3, 3, 0,
		3, 3, 0,
		0, 0, 0
	), 5L, 3L, TRUE, sparse = TRUE)
	expect_is(m, 'dgCMatrix')
	
	expect_identical(duplicated(m), c(FALSE, FALSE, FALSE, TRUE, TRUE))
	expect_identical(duplicated(m, MARGIN = 2), c(FALSE, TRUE, FALSE))
	expect_error(duplicated(m, MARGIN = 3), 'Invalid MARGIN 3')
})

test_that('duplicated.dgCMatrix works with zero columns', {
	m <- sparseMatrix(
		c(1, 1, 1, 2, 1),
		c(1, 2, 4, 4, 6),
		x = c(1, 1, 2, 3, 1))
	expect_is(m, 'dgCMatrix')
	
	expect_identical(which(duplicated(m, MARGIN = 2)), c(2L, 5L, 6L))
})
