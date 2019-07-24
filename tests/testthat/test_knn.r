context('kNN integrity')

test_that('knn works similarly to FNN', {
	skip_if_not_installed('FNN')
	
	data(guo_norm, package = 'destiny')
	e <- t(Biobase::exprs(guo_norm))
	
	r_destiny <- destiny::find_knn(e, 5L)
	r_fnn <- FNN::get.knn(e, 5L)
	
	dimnames(r_destiny$index) <- dimnames(r_destiny$dist) <- NULL
	
	rows_eq <- sapply(seq_len(nrow(e)), function(r) isTRUE(all.equal(r_destiny$index[r, ], r_fnn$nn.index[r, ])))
	expect_lte(sum(!rows_eq), 5L)
	expect_equal(r_destiny$dist[rows_eq, ], r_fnn$nn.dist[rows_eq, ], tolerance = 1e-5)
})

test_that('knnx works similarly to FNN', {
	skip_if_not_installed('FNN')
	
	data(guo_norm, package = 'destiny')
	e <- t(Biobase::exprs(guo_norm))
	nc <- guo_norm$num_cells
	
	r_destiny <- destiny::find_knn(e[nc != 32L,], 5L, query = e[nc == 32L,])
	r_fnn <- FNN::get.knnx(e[nc != 32L,], e[nc == 32L,], 5L)
	
	dimnames(r_destiny$index) <- dimnames(r_destiny$dist) <- NULL
	
	rows_eq <- sapply(seq_len(sum(nc == 32L)), function(r) isTRUE(all.equal(r_destiny$index[r, ], r_fnn$nn.index[r, ])))
	expect_lte(sum(!rows_eq), 30L)
	expect_equal(r_destiny$dist[rows_eq, ], r_fnn$nn.dist[rows_eq, ], tolerance = 1e-5)
})
