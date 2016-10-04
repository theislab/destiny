#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dMatrix ddenseMatrix dsparseMatrix ddiMatrix dgeMatrix dtrMatrix dtpMatrix dsyMatrix dspMatrix dpoMatrix dppMatrix corMatrix Cholesky pCholesky BunchKaufman pBunchKaufman
setClassUnion('characterOrnumericOrNULL', members = c('character', 'numeric', 'NULL'))
setClassUnion('numericOrNULL', members = c('numeric', 'NULL'))
setClassUnion('integerOrNULL', members = c('integer', 'NULL'))
setClassUnion('dMatrixOrNULL', members = c('dMatrix', 'NULL'))

setClassUnion('logicalOrMissing', members = c('logical', 'missing'))
