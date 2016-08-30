bad_default_palette <- c('black', 'red', 'green3', 'blue', 'cyan', 'magenta', 'yellow', 'gray')

#' @importFrom scales brewer_pal
.onLoad <- function(libname, pkgname) {
	if (identical(palette(), bad_default_palette)) {
		better_default_palette <- brewer_pal(palette = 'Set3')
		palette(better_default_palette(12L))
	}
}
