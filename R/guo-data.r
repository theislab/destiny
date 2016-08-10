#' Guo at al. mouse embryonic stem cell qPCR data
#' 
#' Gene expression data of 48 genes and an annotation column \code{$num_cells} containing the cell stage at which the embryos were harvested.
#' 
#' The data is normalized using the mean of two housekeeping genes.
#' The difference between \code{guo} and \code{guo_norm} is the LoD being set to 10 in the former, making it usable with the \code{censor_val} parameter of \link{DiffusionMap}.
#' 
#' @return an \link[Biobase]{ExpressionSet} with 48 features and 428 samples containing qPCR Ct values and a "num.cells" sample annotation.
#' 
#' @aliases data:guo data:guo_norm guo guo_norm
#' @name guo
#' @docType data
#' @author Guoji Guo, Mikael Huss, Guo Qing Tong, Chaoyang Wang, Li Li Sun, Neil D. Clarke, Paul Robson \email{robsonp@@gis.a-star.edu.sg}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S1534580710001103}
#' @keywords data
#' @usage
#' data(guo)
#' data(guo_norm)
#' @format An \link[Biobase]{ExpressionSet} with 48 features, 428 samples and 2 \link[Biobase]{phenoData} annotations.
NULL
