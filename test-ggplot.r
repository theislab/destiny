library(destiny)
library(dplyr)
data(guo_norm)

guo_32   <- guo_norm %>% as('data.frame') %>% filter(num_cells == 32)
dm_no_32 <- guo_norm %>% as('data.frame') %>% filter(num_cells != 32) %>% DiffusionMap

#3d
plot.DiffusionMap(dm_no_32, 1:3, col_by = 'num_cells',                pch = 20)
plot.DiffusionMap(dm_no_32, 1:3, col_by = 'num_cells', ticks =  TRUE, pch = 20)
plot.DiffusionMap(dm_no_32, 1:3, col_by = 'num_cells', axes  = FALSE, pch = 20)
plot.DiffusionMap(dm_no_32, 1:3, col_by = 'num_cells', box   =  TRUE, pch = 20)

#ggplot
plot.DiffusionMap(dm_no_32, 1:2, col_by = 'num_cells')
plot.DiffusionMap(dm_no_32, 1:2, col_by = 'num_cells', ticks =  TRUE)
plot.DiffusionMap(dm_no_32, 1:2, col_by = 'num_cells', axes  = FALSE)
plot.DiffusionMap(dm_no_32, 1:2, col_by = 'num_cells', box   =  TRUE)
