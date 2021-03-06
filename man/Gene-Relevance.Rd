% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene-relevance.r
\docType{class}
\name{GeneRelevance-class}
\alias{GeneRelevance-class}
\alias{gene_relevance}
\alias{gene_relevance,DiffusionMap,missing-method}
\alias{gene_relevance,matrix,dMatrixOrMatrix-method}
\title{Gene relevances for entire data set}
\usage{
gene_relevance(
  coords,
  exprs,
  ...,
  k = 20L,
  dims = 1:2,
  distance = NULL,
  smooth = TRUE,
  remove_outliers = FALSE,
  verbose = FALSE
)

\S4method{gene_relevance}{DiffusionMap,missing}(
  coords,
  exprs,
  ...,
  k = 20L,
  dims = 1:2,
  distance = NULL,
  smooth = TRUE,
  remove_outliers = FALSE,
  verbose = FALSE
)

\S4method{gene_relevance}{matrix,dMatrixOrMatrix}(
  coords,
  exprs,
  ...,
  pcs = NULL,
  knn_params = list(),
  weights = 1,
  k,
  dims,
  distance,
  smooth,
  remove_outliers,
  verbose
)
}
\arguments{
\item{coords}{A \code{\link{DiffusionMap}} object or a cells \eqn{\times} dims \code{\link{matrix}}.}

\item{exprs}{An cells \eqn{\times} genes \code{\link{matrix}}. Only provide if \code{coords} is no \code{\link{DiffusionMap}}.}

\item{...}{Unused. All parameters to the right of the \code{...} have to be specified by name.}

\item{k}{Number of nearest neighbors to use}

\item{dims}{Index into columns of \code{coord}}

\item{distance}{Distance measure to use for the nearest neighbor search.}

\item{smooth}{Smoothing parameters \code{c(window, alpha)} (see \code{\link[smoother]{smth.gaussian}}).
Alternatively \code{\link{TRUE}} to use the \link[smoother]{smoother} \link[smoother:smth.options]{defaults}
or \code{\link{FALSE}} to skip smoothing,}

\item{remove_outliers}{Remove cells that are only within one other cell's nearest neighbor, as they tend to get large norms.}

\item{verbose}{If TRUE, log additional info to the console}

\item{pcs}{A cell \eqn{\times} \code{n_pcs} matrix of principal components to use for the distances.}

\item{knn_params}{A \code{\link{list}} of parameters for \code{\link{find_knn}}.}

\item{weights}{Weights for the partial derivatives. A vector of the same length as \code{dims}.}
}
\value{
A \code{GeneRelevance} object:
}
\description{
The relevance map is cached insided of the \code{\link{DiffusionMap}}.
}
\section{Slots}{

\describe{
\item{\code{coords}}{A cells \eqn{\times} dims \code{\link{matrix}} or \code{\link[Matrix:sparseMatrix-class]{sparseMatrix}}
of coordinates (e.g. diffusion components), reduced to the dimensions passed as \code{dims}}

\item{\code{exprs}}{A cells \eqn{\times} genes matrix of expressions}

\item{\code{partials}}{Array of partial derivatives wrt to considered dimensions in reduced space
(genes \eqn{\times} cells \eqn{\times} dimensions)}

\item{\code{partials_norm}}{Matrix with norm of aforementioned derivatives. (n\_genes \eqn{\times} cells)}

\item{\code{nn_index}}{Matrix of k nearest neighbor indices. (cells \eqn{\times} k)}

\item{\code{dims}}{Column index for plotted dimensions. Can \code{\link{character}}, \code{\link{numeric}} or \code{\link{logical}}}

\item{\code{distance}}{Distance measure used in the nearest neighbor search. See \code{\link{find_knn}}}

\item{\code{smooth_window}}{Smoothing window used (see \code{\link[smoother]{smth.gaussian}})}

\item{\code{smooth_alpha}}{Smoothing kernel width used (see \code{\link[smoother]{smth.gaussian}})}
}}

\examples{
data(guo_norm)
dm <- DiffusionMap(guo_norm)
gr <- gene_relevance(dm)

m <- t(Biobase::exprs(guo_norm))
gr_pca <- gene_relevance(prcomp(m)$x, m)
# now plot them!

}
\seealso{
\link{Gene Relevance methods}, \link{Gene Relevance plotting}: \code{plot_differential_map}/\code{plot_gene_relevance}
}
