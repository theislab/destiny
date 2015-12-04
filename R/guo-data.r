#' Guo at al. mouse embryonic stem cell qPCR data
#' 
#' Gene expression data of 48 genes and an annotation column \code{$num.cells} containing the cell stage at which the embryos were harvested.
#' 
#' The data is normalized using the mean of two housekeeping genes.
#' The difference between \code{guo} and \code{guo.norm} is the LoD being set to 10 in the former, making it usable with the \code{censor.val} parameter of \link{DiffusionMap}.
#' 
#' @return a \eqn{429 \times 49} \code{data.frame} with the first 48 columns being qPCR Ct values and the last one the "divisions" annotation.
#' 
#' @aliases data:guo data:guo.norm guo guo.norm
#' @name guo
#' @docType data
#' @author Guoji Guo, Mikael Huss, Guo Qing Tong, Chaoyang Wang, Li Li Sun, Neil D. Clarke, Paul Robson \email{robsonp@@gis.a-star.edu.sg}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S1534580710001103}
#' @keywords data
#' @usage
#' data(guo)
#' data(guo.norm)
#' @format A data frame with 492 rows and 49 variables, the first 48 of which are Ct values for gene expression data
NULL
