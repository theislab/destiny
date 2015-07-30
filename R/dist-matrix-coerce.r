setOldClass('dist')

as.Matrix.dist <- function(from) {
	s <- attr(from, 'Size')
	i <- rev(abs(sequence(seq.int(s - 1L)) - s) + 1L)
	j <- rep.int(seq_len(s - 1L), rev(seq_len(s - 1L)))
	sparseMatrix(i, j, x = unclass(from), dims = c(s, s), symmetric = TRUE)
}

setAs('dist', 'Matrix',        as.Matrix.dist)
setAs('dist', 'sparseMatrix',  as.Matrix.dist)
setAs('dist', 'CsparseMatrix', as.Matrix.dist)
setAs('dist', 'dsCMatrix',     as.Matrix.dist)
