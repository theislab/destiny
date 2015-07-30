getVigSources <- function(vigdir) {
	pattern = '[.]Rmd$|[.]Rnw$|[.]Rrst$|[.]Rhtml$|[.]Rtex$'
	builder.pkg <- read.dcf(file.path(vigdir, '..', 'DESCRIPTION'), 'VignetteBuilder')[[1]]
	if (!is.na(builder.pkg)) {
		builders <- tryCatch(tools::vignetteEngine(package = builder.pkg), error = function(e) NULL)
		if (!is.null(builders)) {
			builders.pattern <- paste(lapply(builders, function(b) b$pattern), collapse = '|')
			pattern <- paste(pattern, builders.pattern, sep = '|')
		}
	}
	dir(vigdir, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
}