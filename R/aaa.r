bad_default_palette <- c('black', 'red', 'green3', 'blue', 'cyan', 'magenta', 'yellow', 'gray')
#this is a reshuffled Set3, but RColorBrewer doesn’t work in .onLoad due to using rgb() unconditionally
gud_default_palette <- c('#8DD3C7', '#FFED6F', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#BC80BD', '#FCCDE5', '#D9D9D9', '#CCEBC5', '#FFFFB3')

# do not use .Call.graphics like grDevices::palette
set_palette <- function(v) .Call(grDevices:::C_palette, v)

.onLoad <- function(libname, pkgname) {
	if (identical(palette(), bad_default_palette)) {
		set_palette(gud_default_palette)
	}
}
