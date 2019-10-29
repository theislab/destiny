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
	if (FALSE) {
		library(ggplot2)
		library(tidyr)
		library(dplyr)
		gg_heatmap <- cbind(
			as.data.frame(scales::rescale(gr@exprs)),
			as.data.frame(scales::rescale(gr@partials_norm)) %>% rename_all(~ paste(., 'norm')),
			as.data.frame(gr@coords) %>% select(DC1, DC2)
		) %>% mutate(Cell = seq_len(n()), OrderDC = order(DC1 + DC2)) %>%
			gather(Gene, Expr, -Cell, -DC1, -DC2, -OrderDC) %>%
			ggplot(aes(Gene, OrderDC, fill = Expr)) +
				geom_tile() +
				scale_x_discrete(expand = c(0, 0), position = 'top') +
				scale_y_reverse(expand = c(0, 0)) +
				scale_fill_viridis_c() +
				theme(axis.ticks.x = element_blank())
		print(gg_heatmap)
		print(plot_differential_map(gr, gene = colnames(test_data)[1:4]))
		gr_plot %>% with(data.frame(ids, scores)) %>% print()
		print(plot(dm, 1:2, col_by = 'Type'))
		print(gr_plot)
	}
	skip('Not yet stable')
	expect_identical(gr_plot$ids[1:2], c('D', 'C'))
})
