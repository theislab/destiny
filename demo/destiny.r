library(destiny)
library(RColorBrewer)

data(guo)

palette(brewer.pal(8L, 'Dark2'))




sigmas <- find.sigmas(guo, censor.val = 10, censor.range = c(10, 40),
                      verbose = FALSE)
par(lwd = 3)
plot(sigmas,
		 col           = palette()[[1]],
		 highlight.col = palette()[[4]],
		 line.col      = palette()[[6]],
     rect.options = list(border = 'white'))




diff.guo <- DiffusionMap(guo, sigmas, censor.val = 10, censor.range = c(10, 40),
                         verbose = FALSE)
plot(diff.guo,
     col = guo$num.cells, pch = 20)

#library(rgl)
#plot3d(eigenvectors(diff.guo)[, 1:3], col = guo$num.cells)
