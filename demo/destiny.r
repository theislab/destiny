library(destiny)

data(guo)

dark2 <- scales::brewer_pal(palette = 'Dark2')
palette(dark2(8L))




dm_guo <- DiffusionMap(guo, verbose = FALSE,
                       censor_val = 10, censor_range = c(10, 40))
plot(dm_guo,
     col = guo$num_cells, pch = 20)




sigmas <- find_sigmas(guo, verbose = FALSE,
                      censor_val = 10, censor_range = c(10, 40))
par(lwd = 3)
plot(sigmas,
    col           = palette()[[1]],
    col_highlight = palette()[[4]],
    col_line      = palette()[[6]])




dm_guo_global <- DiffusionMap(guo, sigmas, verbose = FALSE,
                              censor_val = 10, censor_range = c(10, 40))
plot(dm_guo_global,
     col = guo$num_cells, pch = 20)
