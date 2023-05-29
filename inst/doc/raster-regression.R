## ----setup, message = FALSE, warning = FALSE, eval = TRUE---------------------
library(geostan)
library(sf)

## ----eval = TRUE--------------------------------------------------------------
row <- 40
col <- 30
c(N <- row * col)
sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(col,0), c(col,row), c(0,0)))))
grid <- st_make_grid(sfc, cellsize = 1, square = TRUE)
grid <- st_as_sf(grid)
W <- shape2mat(grid, style = "W", queen = FALSE)
grid$z <- sim_sar(w = W, rho = 0.9)
grid$y <- -0.5 * grid$z + sim_sar(w = W, rho = .9, sigma = .3)

## ----fig.width = 3.25, fig.height = 3, fig.align = 'center'-------------------
plot(grid[,'z'])

## ----eval = FALSE-------------------------------------------------------------
#  fit <- stan_sar(y ~ z, data = grid, C = W)

## ----eval = FALSE-------------------------------------------------------------
#  C <- shape2mat(grid, style = "B", queen = FALSE)
#  car_list <- prep_car_data(C, "WCAR")
#  fit <- stan_car(y ~ z, data = grid, car_parts = car_list)

## ----eval = TRUE--------------------------------------------------------------
sar_list <- prep_sar_data2(row = row, col = col)
fit <- stan_sar(
                y ~ z, 
                data = grid,
		centerx = TRUE,
		sar_parts = sar_list,
		iter = 500,
		chains = 4,
		slim = TRUE #,
		# cores = 4, # for multi-core processing
		)
print(fit)		

