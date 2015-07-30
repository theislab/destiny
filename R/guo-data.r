#' Guo at al. mouse embryonic stem cell qPCR data
#' 
#' Gene expression data of 48 genes and an annotation column \code{$division} containing the number of divisions after which the cells were harvested.
#' 
#' The data is normalized using the mean of two housekeeping genes and the LoD is set to 15.
#' 
#' @return a \eqn{429 \times 49} \code{data.frame} with the first 48 columns being qPCR Ct values and the last one the "divisions" annotation.
#' 
#' @aliases data(guo) data-guo guo
#' @name data-guo
#' @docType data
#' @author Guoji Guo, Mikael Huss, Guo Qing Tong, Chaoyang Wang, Li Li Sun, Neil D. Clarke, Paul Robson \email{robsonp@@gis.a-star.edu.sg}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S1534580710001103}
#' @keywords data
#' @usage data(guo)
#' @format A data frame with 492 rows and 49 variables, the first 48 of which are Ct values for gene expression data
NULL
