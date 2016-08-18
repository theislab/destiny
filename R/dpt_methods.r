#' @export
as.data.frame.DPT <- function(x, row.names = NULL, optional = FALSE, ...) data.frame(
	Branch = x@branch,
	Parent = x@parent,
	Tip    = lWhich(x@tips, len = length(x@branch)),
	DPT    = x@dpt[, 1],  #DPT
	DPT    = x@dpt[,-1])  #DPT.1, DPT.2, ...