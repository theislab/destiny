# Smooth Using Gaussian Window
# 
# The specific function for smoothing using the gaussian window function
# @param x numeric vector of values to smooth, error will be thrown if not provided.
# @param window the length of the smoothing window, if an integer, represents
# number of items, else, if a value between \code{0} and \code{1}, represents the 
# proportion of the input vector
# @param alpha parameter to determine the breadth of the gaussian window, yielding more or less sensitive 
# smoothing characteristics
# @param ... not used
# @param tails Logical value as to whether the tail regions should be included or not.
smth.gaussian <- function(
	x       = stop("Numeric Vector 'x' is Required"),
  window  = getOption('smoother.window', 0.1),
  alpha   = getOption('smoother.gaussianwindow.alpha', 2.5),
  ...,
  tails   = getOption('smoother.tails', FALSE)){
	
	if (!is.numeric(x) | !is.numeric(alpha)) {
		stop("argument 'x' and 'alpha' must be numeric", call. = FALSE)
	}
	windowLength = .determineWindowLength(x, window)
	makeWindow = function(w, a) {
		hw = abs(w/2)
		e = exp(1)
		a = abs(a)
		ret = sapply(c(0:(w - 1)), function(x) {
			n = x - as.integer(hw)
			k = -0.5 * (a * n/hw)^2
			e^k
		})
		ret
	}
	w = makeWindow(windowLength, alpha[1])
	sizeW = length(w)
	sizeD = length(x)
	w = .normalize(w)
	hkwL = as.integer(sizeW/2)
	hkwR = sizeW - hkwL
	ret = sapply(c(1:sizeD), function(i) {
		ix.d = c((i - hkwL):(i + hkwR - 1))
		ix.w = which(ix.d %in% 1:sizeD)
		ix.d = ix.d[ix.w]
		W.nm = if(length(ix.w) != sizeW){.normalize(w[ix.w])}else{w}
		D.nm = x[ix.d]
		as.numeric(D.nm %*% W.nm)
	})
	if (!tails) {
		ret[c(1:hkwL, (sizeD - hkwR + 1):sizeD)] = NA
	}
	ret
}

.determineWindowLength = function (x, w) {
	if (!is.numeric(x) | !is.numeric(w)) {
		stop("Arguments 'x' and 'w' must be numeric")
	}
	l = length(x)
	w.orig = w
	if (w > 0 && w < 1) {
		w = w * l
	}
	ret = as.integer(max(abs(w[1]), 1))
	if (length(x) <= ret | ret < 1) {
		stop(paste("Resultant window length is out of range (", 
							 ret, "), must be >= 1", sep = ""), call. = FALSE)
	}
	return(ret)
}

.normalize = function(x) {
	if (!is.numeric(x)) {
		stop("argument 'x' must be numeric")
	}
	if (length(x) == 1) {
		return(1)
	}
	total = sum(x)
	if (total == 0) {
		stop("argument 'x' must not sum to zero")
	}
	return(x/total)
}

