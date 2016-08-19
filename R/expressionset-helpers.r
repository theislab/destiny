#' @include utils.r
NULL

#' Convert object to \link[Biobase]{ExpressionSet} or read it from a file
#' 
#' These functions present quick waya to create \link[Biobase]{ExpressionSet} objects.
#' 
#' They work by using all continuous (double) columns as expression data, and all others as sample annotations.
#' 
#' @param x  \link[base]{data.frame} to convert to an \link[Biobase]{ExpressionSet}.
#' 
#' @return an \link[Biobase]{ExpressionSet} object
#' 
#' @examples
#' library(Biobase)
#' df <- data.frame(Time  = seq_len(3), #integer column
#'                  Actb  = c(0.05, 0.3, 0.8),
#'                  Gapdh = c(0.2, 0.03, 0.1))
#' set <- as.ExpressionSet(df)
#' rownames(exprs(set)) == c('Actb', 'Gapdh')
#' phenoData(set)$Time == 1:3
#' 
#' @seealso \link[utils]{read.table} on which \code{read.ExpressionSet} is based, and \link[Biobase]{ExpressionSet}.
#' 
#' @aliases as.ExpressionSet as.ExpressionSet-method as.ExpressionSet,data.frame-method read.ExpressionSet
#' 
#' @importFrom methods setGeneric .valueClassTest
#' @name ExpressionSet helpers
#' @export
setGeneric('as.ExpressionSet', function(x, ...) standardGeneric('as.ExpressionSet'), valueClass = 'ExpressionSet')

#' @param annotation_cols  The data.frame columns used as annotations. All others are used as expressions. (Logical, character or numerical index array)
#' 
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame phenoData
#' @name ExpressionSet helpers
#' @export
setMethod('as.ExpressionSet', 'data.frame', function(x, annotation_cols = !sapply(x, is.double)) {
	if (!is.logical(annotation_cols))
		annotation_cols <- lWhich(annotation_cols, names(x))
	
	assayData <- t(as.matrix(x[!annotation_cols]))
	phenoData <- AnnotatedDataFrame(x[annotation_cols])
	
	ExpressionSet(assayData, phenoData)
})

#' @param file    File path to read ASCII data from
#' @param header  Specifies if the file has a header row.
#' @param ...     Additional parameters to \link[utils]{read.table}
#' 
#' @importFrom utils read.table
#' @name ExpressionSet helpers
#' @export
read.ExpressionSet <- function(file, header = TRUE, ...) {
	as.ExpressionSet(read.table(file, header, ...))
}
