context('kNN integrity')

test_that('knn works the same as in FNN', {
	skip_if_not_installed('FNN')
	
	data(guo_norm, package = 'destiny')
	e <- t(Biobase::exprs(guo_norm))
	
	r_destiny <- destiny::find_knn(e, 5L)
	r_fnn <- FNN::get.knn(e, 5L)
	
	dimnames(r_destiny$index) <- dimnames(r_destiny$dist) <- NULL
	expect_equal(r_destiny$index, r_fnn$nn.index)
	expect_equal(r_destiny$dist,  r_fnn$nn.dist)
})

test_that('knnx works the same as in FNN', {
	skip_if_not_installed('FNN')
	
	data(guo_norm, package = 'destiny')
	e <- t(Biobase::exprs(guo_norm))
	nc <- guo_norm$num_cells
	
	r_destiny <- destiny::find_knn(e[nc != 32L,], 5L, query = e[nc == 32L,])
	r_fnn <- FNN::get.knnx(e[nc != 32L,], e[nc == 32L,], 5L)
	
	dimnames(r_destiny$index) <- dimnames(r_destiny$dist) <- NULL
	expect_equal(r_destiny$index, r_fnn$nn.index)
	expect_equal(r_destiny$dist,  r_fnn$nn.dist)
})
