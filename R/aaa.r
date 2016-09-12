bad_default_palette <- c('black', 'red', 'green3', 'blue', 'cyan', 'magenta', 'yellow', 'gray')
#this is Set3, but RColorBrewer doesn’t work in .onLoad due to using rgb() unconditionally
gud_default_palette <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F')

.onLoad <- function(libname, pkgname) {
	if (identical(palette(), bad_default_palette)) {
		palette(gud_default_palette)
	}
}
