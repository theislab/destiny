context('utils')

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
	expect_identical(duplicated(t(m)), c(FALSE, TRUE, FALSE))
})
