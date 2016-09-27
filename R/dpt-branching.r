auto_branch <- function(acc_trans, stats, w_width, nmin = 10L, gmin = 1.1) {
	n <- ncol(acc_trans)
	
	stopifnot(n >= nmin)
	stopifnot(stats$g >= gmin)
	
	# initialize one level (branch, tips) and three branches (dpt)
	dpt <- dpt_from_tips(acc_trans, stats$tips)
	branches <- cut_branches(dpt, w_width)  # list ov vectors of numeric indices
	branch <- matrix(idx_list_to_vec(branches, n), n, 1L)
	tips <- matrix(logical(n), n, 1L)
	tips[stats$tips, 1L] <- TRUE
	
	subs <- lapply(seq_len(3L), function(s) {
		idx_sub <- branches[[s]]
		i <- stats$tips[[s]]
		if (length(idx_sub) < nmin || !i %in% idx_sub)  # Tip cells can end up undecided!
			return(NULL)
		
		sub_props <- as(acc_trans[idx_sub, idx_sub, drop = FALSE], 'symmetricMatrix')
		sub_root  <- abs_idx_to_sub_idx(idx_sub, i)
		sub_stats <- tipstats(sub_props, sub_root)
		if (sub_stats$g < gmin)
			return(NULL)
		
		auto_branch(sub_props, sub_stats, w_width, nmin, gmin)
	})
	
	# add dpt columns to dpt and a level column to branch/tips if any branch was subdivided
	nonnull_subs <- vapply(subs, Negate(is.null), logical(1L))
	if (any(nonnull_subs)) {
		n_sublevels <- do.call(max, lapply(subs[nonnull_subs], function(s) ncol(s$branch)))
		
		branch <- cbind(branch, matrix(NA_integer_, n, n_sublevels))
		tips   <- cbind(tips,   matrix(NA,          n, n_sublevels))
		
		for (s in which(nonnull_subs)) {
			sub <- subs[[s]]
			idx_sub <- branches[[s]]
			
			d_new <- matrix(NA_real_, n, ncol(sub$dpt))
			d_new[idx_sub, ] <- sub$dpt
			dpt <- cbind(dpt, d_new)
			
			idx_newcol <- seq.int(ncol(branch) - n_sublevels + 1L, length.out = ncol(sub$branch))
			stopifnot(ncol(sub$branch) == ncol(sub$tips))
			
			branch_offset <- max(branch, na.rm = TRUE)
			branch[idx_sub, idx_newcol] <- sub$branch + branch_offset
			tips[  idx_sub, idx_newcol] <- sub$tips
		}
	}
	
	stopifnot(ncol(branch) == ncol(tips))
	list(dpt = dpt, branch = branch, tips = tips)
}


idx_list_to_vec <- function(idx_list, n) {
	v <- rep(NA_integer_, n)
	for (i in seq_along(idx_list))
		v[idx_list[[i]]] <- i
	v
}


# v <- c(F,T,F,T,T,F)
# abs_idx_to_sub_idx(which(v), c(2, 5)) -> c(1, 3)
abs_idx_to_sub_idx <- function(idx_abs, i) {
	n_new <- length(idx_abs)
	idx_sub <- rep(NA_integer_, max(i, na.rm = TRUE))
	idx_sub[idx_abs] <- seq_len(n_new)
	idx_sub[i]
}


dpt_from_tips <- function(acc_trans, tips) {
	n <- ncol(acc_trans)
	vapply(tips, function(cell) dpt_to_cell(acc_trans, cell), double(n))
}


cut_branches <- function(dpt, w_width) {
	bid <- apply(dpt, 2, order)
	branches <- lapply(seq_len(3), function(b) branchcut(dpt, bid, b, w_width))
	organize_branches(branches)
}


#' @importFrom smoother smth.gaussian
branchcut <- function(dpt, bid, b, w_width) {
	n <- nrow(bid)
	all_branches <- seq_len(3L)
	
	# sanity checks
	stopifnot(b %in% all_branches)
	stopifnot(ncol(dpt) == 3L, ncol(bid) == 3L)
	stopifnot(nrow(dpt) == n)
	stopifnot(is.double(dpt), is.integer(bid))
	
	# find cell indexes per branch 
	other <- all_branches[all_branches != b]
	b1 <- other[[1L]]
	b2 <- other[[2L]]
	
	# DPT for other branches, sorted by b3
	b3_idxs <- bid[, b]
	dpt1 <- dpt[b3_idxs, b1]
	dpt2 <- dpt[b3_idxs, b2]
	
	kcor <- vapply(seq_len(n - 1L), function(s1) {
		s2 <- s1 + 1L
		l <- seq_len(s1)
		r <- seq(s2, n)
		
		k_l <- kendall_finite_cor(dpt1[l], dpt2[l], dpt1[[s2]], dpt2[[s2]])
		k_r <- kendall_finite_cor(dpt1[r], dpt2[r], dpt1[[s1]], dpt2[[s1]])
		
		k_l/s1 - k_r/(n - s1)
	}, double(1L))
	
	kcor <- smth.gaussian(kcor, w_width)
	cut <- which.max(kcor)
	
	b3_idxs[seq_len(cut)]
}


# This function organizes the cell arrays branches[[i]]
# to build Branch labels (length <- number of cells) indicating the branch each cell belongs to.
# Cells which are assigned to more than one branch in dpt as well
# as cells which are not assigned to any branch are defined as undeceided (label NA)
#' @importFrom utils combn
organize_branches <- function(branches) {
	n <- do.call(max, branches)
	
	intersect_branches <- function(bs) intersect(branches[[bs[[1L]]]], branches[[bs[[2L]]]])
	branch_intersections <- lapply(combn(3L, 2L, simplify = FALSE), intersect_branches)
	inters <- Reduce(union, branch_intersections, integer())
	
	lapply(branches, function(b) setdiff(b, inters))
}


# given two orderings b1 and b2, compute the delta (e.i. finite) Kendall
# correlation of adding a new cell with bnew1, bnew2 to the orderings.
kendall_finite_cor <- function(b1, b2, new1, new2) {
	b11 <- numeric(length(b1))
	b11[b1 >= new1] <- 1
	b11[b1 <  new1] <- -1
	
	b22 <- numeric(length(b2))
	b22[b2 >= new2] <- 1
	b22[b2 <  new2] <- -1
	
	as.double(b11 %*% b22)  # strip dims
}
