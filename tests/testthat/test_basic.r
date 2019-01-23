context('basic functionality')

test_that('DiffusionMap works', {
	data(guo_norm)
	dm <- DiffusionMap(guo_norm)
	expect_identical(dm@distance, 'euclidean')
})

test_that('Eigenvectors retain cell names', {
	data(guo_norm)
	dm <- DiffusionMap(guo_norm)
	expect_identical(rownames(dm@eigenvectors), Biobase::sampleNames(guo_norm))
})

test_that('DM works with dist matrices', {
	data(guo_norm)
	dists <- as(dist(t(Biobase::exprs(guo_norm))), 'sparseMatrix')
	rownames(dists) <- colnames(dists) <- Biobase::sampleNames(guo_norm)
	dm <- DiffusionMap(distance = dists)
	expect_identical(rownames(dm@eigenvectors), rownames(dists))
})
