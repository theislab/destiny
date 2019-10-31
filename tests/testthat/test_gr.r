context('gene relevance')

smps <- 120L
part <- 1/4
test_data <- local({
	d <- data.frame(
	  A = c(                   seq(1,  0, length = smps * part), rep(0, smps*3*part)),
		B = c(                   seq(0,  1, length = smps*2*part), rep(1, smps*2*part)),
		C = c(rep(0, smps*part), seq(0,  1, length = smps*2*part), rep(0, smps * part)),
		D = c(rep(0, smps*part), seq(0, .5, length = smps * part), rep(0, smps * part), seq(.5, 1, length = smps*part)), #  + rnorm(smps, 0, .05)
		stringsAsFactors = FALSE
	)
	d$Cell <- seq_len(nrow(d))
	d$Type <- as.integer(ceiling(d$Cell / (smps/4)))
	d
})

test_that('Gene Relevance returns the right genes', {
	set.seed(0)
	dm <- DiffusionMap(test_data, distance = 'rankcor')
	gr <- gene_relevance(dm,  smooth = FALSE)
	gr_plot <- plot_gene_relevance(gr, iter_smooth = 0, n_top = 1)
	skip('Not yet stable')
	expect_identical(gr_plot$ids[1:2], c('D', 'C'))
})
