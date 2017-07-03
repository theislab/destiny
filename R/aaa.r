bad_default_palette <- c('black', 'red', 'green3', 'blue', 'cyan', 'magenta', 'yellow', 'gray')
#this is a reshuffled Set3, but RColorBrewer doesn’t work in .onLoad due to using rgb() unconditionally
gud_default_palette <- c('#8DD3C7', '#FFED6F', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#BC80BD', '#FCCDE5', '#D9D9D9', '#CCEBC5', '#FFFFB3')

# do not use .Call.graphics like grDevices::palette
#' @importFrom utils getFromNamespace
set_palette <- function(v) {
	cp <- getFromNamespace('C_palette', 'grDevices')
	do.call(.Call, list(cp, v))
}

.onLoad <- function(libname, pkgname) {
	if (identical(palette(), bad_default_palette)) {
		set_palette(gud_default_palette)
	}
}

# work around scatterplot3d bug. use function from this package so `linepad` isn’t unset
# #' @importFrom scatterplot3d scatterplot3d
#' @importFrom grDevices xyz.coords
#' @importFrom graphics lines mtext strheight strwidth title
#' @importFrom stats coef
linepad <- -.5
scatterplot3d <- scatterplot3d::scatterplot3d
environment(scatterplot3d) <- environment(set_palette)
