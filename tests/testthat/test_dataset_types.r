library(Biobase)

test_nobss <- 4L
test_nfeat <- 3L

test_obss <- paste0('c', seq_len(test_nobss))
test_feat <- paste0('g', seq_len(test_nfeat))

test_matrix <- matrix(runif(test_nobss * test_nfeat), test_nobss, test_nfeat, dimnames = list(test_obss, test_feat))
test_df     <- data.frame(test_matrix, cm1 = LETTERS[seq_len(test_nobss)])
test_es     <- ExpressionSet(
	t(test_matrix),
	AnnotatedDataFrame(data.frame(cm1 = LETTERS[seq_len(test_nobss)], row.names = test_obss)),
	AnnotatedDataFrame(data.frame(gm1 = letters[seq_len(test_nfeat)], row.names = test_feat)))


test_that('The helpers work with matrix data', {
	expect_identical(dataset_extract_doublematrix(test_matrix), test_matrix)
	expect_identical(dataset_n_observations      (test_matrix), test_nobss)
	expect_identical(dataset_n_features          (test_matrix), test_nfeat)
	expect_identical(dataset_to_df               (test_matrix), test_df[, 1:3])
	expect_identical(dataset_names               (test_matrix), test_feat)
	expect_identical(dataset_get_feature   (test_matrix, 'g2'), test_matrix[, 'g2'])
})


test_that('The helpers work with data.frame data', {
	expect_identical(dataset_extract_doublematrix(test_df), test_matrix)
	expect_identical(dataset_n_observations      (test_df), test_nobss)
	expect_identical(dataset_n_features          (test_df), test_nfeat)
	expect_identical(dataset_to_df               (test_df), test_df)
	expect_identical(dataset_names               (test_df), c(test_feat, 'cm1'))
	expect_identical(dataset_get_feature   (test_df, 'g3'), test_df$g3)
	expect_identical(dataset_get_feature   (test_df, 'cm1'), test_df$cm1)
})


test_that('The helpers work with ExpressionSet data', {
	expect_identical(dataset_extract_doublematrix(test_es), test_matrix)
	expect_identical(dataset_n_observations      (test_es), test_nobss)
	expect_identical(dataset_n_features          (test_es), test_nfeat)
	expect_identical(dataset_to_df               (test_es), test_df)
	expect_identical(dataset_names               (test_es), c(test_feat, 'cm1'))
	expect_identical(dataset_get_feature   (test_es, 'g1'), exprs(test_es)['g1', ])
	expect_identical(dataset_get_feature   (test_es, 'cm1'), test_es$cm1)
})
