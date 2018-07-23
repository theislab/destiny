context('ggplot')

data(guo_norm)

get_geom <- function(p, name) {
	
}

guo_df <- as(guo_norm, 'data.frame')
guo_32    <- guo_df[guo_df$num_cells == 32, ]
guo_no_32 <- guo_df[guo_df$num_cells != 32, ]
dm_no_32 <- DiffusionMap(guo_no_32)

#3d. no idea how to test this
if (FALSE) {
	plot.DiffusionMap(dm_no_32, 1:3, col_by = 'num_cells',                pch = 20)
	plot.DiffusionMap(dm_no_32, 1:3, col_by = 'num_cells', ticks =  TRUE, pch = 20)
	plot.DiffusionMap(dm_no_32, 1:3, col_by = 'num_cells', axes  = FALSE, pch = 20)
	plot.DiffusionMap(dm_no_32, 1:3, col_by = 'num_cells', box   =  TRUE, pch = 20)
}

test_that('ggplot plots have the ticks/boxes they should have', {
	p1 <- plot.DiffusionMap(dm_no_32, 1:2, col_by = 'num_cells')
	p2 <- plot.DiffusionMap(dm_no_32, 1:2, col_by = 'num_cells', ticks =  TRUE)
	p3 <- plot.DiffusionMap(dm_no_32, 1:2, col_by = 'num_cells', axes  = FALSE)
	p4 <- plot.DiffusionMap(dm_no_32, 1:2, col_by = 'num_cells', box   =  TRUE)
	
	# check range_frame
	expect_identical(length(p1$layers), 2L)
	expect_identical(length(p2$layers), 2L)
	expect_identical(length(p3$layers), 1L)
	expect_identical(length(p4$layers), 2L)
	
	# check ticks
	expect_identical(class(p1$theme$axis.ticks)[[1L]], 'element_blank')
	expect_identical(class(p2$theme$axis.ticks)[[1L]], 'element_line')
	expect_identical(class(p3$theme$axis.ticks)[[1L]], 'element_blank')
	expect_identical(class(p4$theme$axis.ticks)[[1L]], 'element_blank')
	
	# check box
	expect_identical(class(p1$theme$panel.border)[[1L]], 'element_blank')
	expect_identical(class(p2$theme$panel.border)[[1L]], 'element_blank')
	expect_identical(class(p3$theme$panel.border)[[1L]], 'element_blank')
	expect_identical(class(p4$theme$panel.border)[[1L]], 'element_rect')
})
