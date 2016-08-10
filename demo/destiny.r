library(destiny)
library(RColorBrewer)

data(guo)

palette(brewer.pal(8L, 'Dark2'))




dm_guo <- DiffusionMap(guo, censor_val = 10, censor_range = c(10, 40),
											 verbose = FALSE)
plot(dm_guo,
		 col = guo$num_cells, pch = 20)




sigmas <- find_sigmas(guo, censor_val = 10, censor_range = c(10, 40),
                      verbose = FALSE)
par(lwd = 3)
plot(sigmas,
		 col           = palette()[[1]],
		 highlight_col = palette()[[4]],
		 line_col      = palette()[[6]],
     rect_options = list(border = 'white'))




dm_guo_global <- DiffusionMap(guo, sigmas, censor_val = 10, censor_range = c(10, 40),
                       verbose = FALSE)
plot(dm_guo_global,
     col = guo$num_cells, pch = 20)




#library(rgl)
#plot3d(eigenvectors(dm_guo)[, 1:3], col = guo$num_cells)
