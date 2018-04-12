test_that('DiffusionMap works', {
	data(guo_norm)
	dm <- DiffusionMap(guo_norm)
	expect_identical(dm@distance, 'euclidean')
})
