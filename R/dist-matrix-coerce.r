#'@importFrom methods setOldClass
setOldClass('dist')

#' @importFrom Matrix sparseMatrix
as.Matrix.dist <- function(from) {
	s <- attr(from, 'Size')
	i <- rep.int(seq_len(s - 1L), rev(seq_len(s - 1L)))
	j <- rev(abs(sequence(seq.int(s - 1L)) - s) + 1L)
	sparseMatrix(i, j, x = unclass(from), dims = c(s, s), symmetric = TRUE)
}

#'@importFrom methods setAs
setAs('dist', 'Matrix',          as.Matrix.dist)
setAs('dist', 'sparseMatrix',    as.Matrix.dist)
setAs('dist', 'CsparseMatrix',   as.Matrix.dist)
setAs('dist', 'dsCMatrix',       as.Matrix.dist)
setAs('dist', 'symmetricMatrix', as.Matrix.dist)
