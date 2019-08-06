context('PCA')

library(Biobase)
library(SingleCellExperiment)

test_nobss <- 30L
test_nfeat <- 20L
test_n_pcs <- 4L

test_obss <- paste0('c', seq_len(test_nobss))
test_feat <- paste0('g', seq_len(test_nfeat))

test_matrix <- matrix(runif(test_nobss * test_nfeat), test_nobss, test_nfeat, dimnames = list(test_obss, test_feat))
test_matrix[test_matrix < mean(test_matrix)] <- 0
test_matrix_sparse <- as(test_matrix, 'sparseMatrix')

test_se <- SingleCellExperiment(
	assays = list(logcounts = t(test_matrix)),
	colData = DataFrame(cm1 = LETTERS[seq_len(test_nobss)]),
	rowData = DataFrame(gm1 = letters[seq_len(test_nfeat)]))
test_se_sparse <- test_se
assay(test_se_sparse) <- t(test_matrix_sparse)


test_that('PCA works sparse and dense data', {
	# Each PC can be flipped individually,
	# so the correct alignment would be to see which ones need to be flipped
	# Using abs() is simpler for comparison, but destroys the real PCs.
	pcs_dense  <- abs(pca_scores(test_matrix, test_n_pcs))
	pcs_sparse <- abs(pca_scores(test_matrix_sparse, test_n_pcs))
	expect_equal(pcs_dense, pcs_sparse, tolerance = 1e-5)
})


test_that('PCA with sparse data does not densify', with_mock(
	prcomp   = function(...) stop('prcomp should not be called'),
	princomp = function(...) stop('princomp should not be called'),
	as.matrix = function(...) stop('as.matrix should not be called'),
	{
		pca_scores(test_matrix_sparse, test_n_pcs)
		succeed('The matrix was not densified')
	}
))
