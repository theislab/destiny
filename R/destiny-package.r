#' Create and plot diffusion maps
#' 
#' The main function is \code{\link{DiffusionMap}}, which returns an object you can \code{\link{plot}} (\code{\link{plot.DiffusionMap}} is then called).
#' 
#' The \code{sigma} parameter for \code{\link{DiffusionMap}} can be determined using \code{\link{find.sigmas}} if your dataset is small enough.
#' 
#' @examples
#' demo(destiny, ask = FALSE)
#' 
#' @docType package
#' @aliases destiny-package
#' @name destiny
#' 
## Make sure Rcpp and RcppEigen are loaded, and the checks don't complain about methods
#' @importFrom Rcpp evalCpp
#' @importFrom RcppEigen RcppEigen.package.skeleton
#' @importFrom methods .valueClassTest validObject new is as
NULL
