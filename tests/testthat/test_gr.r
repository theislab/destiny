smps <- 120L
part <- 1/4
test_data <- local({
	d <- data.frame(
	 #A = c(                   seq(1,  0, length = smps * part), rep(0, smps*3*part)),
		B = c(                   seq(0,  1, length = smps*2*part), rep(1, smps*2*part)),
		C = c(rep(0, smps*part), seq(0,  1, length = smps*2*part), rep(0, smps * part)),
		D = c(rep(0, smps*part), seq(0, .5, length = smps * part), rep(0, smps * part), seq(.5, 1, length = smps*part)), #  + rnorm(smps, 0, .05)
		stringsAsFactors = FALSE
	)
	d$Cell <- seq_len(nrow(d))
	d$Type <- as.integer(ceiling(d$Cell / (smps/4)))
	d
})

if (FALSE) {  # Interactive test
	library(tidyr)
	test_data %>%
		gather(Gene, Expr, -Cell, -Type) %>%
		ggplot(aes(Gene, Cell, fill = Expr)) +
			geom_tile() +
			scale_x_discrete(expand = c(0, 0), position = 'top') +
			scale_y_reverse(expand = c(0, 0)) +
			theme(axis.ticks.x = element_blank())
}

test_that('Gene Relevance returns the right genes', {
	dm <- DiffusionMap(test_data)
	gr <- gene_relevance(dm)
	gr_plot <- plot(gr)
	if (FALSE) {
		print(plot(dm, 1:2, col_by = 'Type'))
		print(gr_plot)
	}
	expect_identical(gr_plot$ids, c('C', 'D'))
})
