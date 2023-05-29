## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      eval = TRUE, 
                      fig.align = "center",
                      fig.width = 3.5,
                      fig.height = 3
                      )

## ----message = FALSE, warning = FALSE-----------------------------------------
library(geostan)
data(georgia)

## -----------------------------------------------------------------------------
georgia$insurance <- georgia$insurance / 100
georgia$insurance.se <- georgia$insurance.se / 100

## -----------------------------------------------------------------------------
max(georgia$insurance + georgia$insurance.se*2)
min(georgia$insurance - georgia$insurance.se*2)

## -----------------------------------------------------------------------------
SE <- data.frame(insurance = georgia$insurance.se)

## -----------------------------------------------------------------------------
C <- shape2mat(georgia, "B")
cars <- prep_car_data(C)

## -----------------------------------------------------------------------------
ME_list <- prep_me_data(se = SE,
                        car_parts = cars,
                        logit = TRUE
                        )

## ----eval = FALSE, message = FALSE--------------------------------------------
#  fit <- stan_car(deaths.male ~ offset(log(pop.at.risk.male)) + insurance,
#                  centerx = TRUE,
#                  data = georgia,
#                  family = poisson(),
#                  car_parts = cars)

## ----eval = TRUE, message = FALSE---------------------------------------------
fit <- stan_car(deaths.male ~ offset(log(pop.at.risk.male)) + insurance,
                centerx = TRUE,
                ME = ME_list,
                family = poisson(),
                data = georgia, 
                car_parts = cars,
                iter = 650, # for demo speed 
                refresh = 0, # minimizes printing
                )

print(fit)

## ----eval = FALSE-------------------------------------------------------------
#  fit <- stan_car(deaths.male ~ offset(log(pop.at.risk.male)) + insurance,
#                  centerx = TRUE,
#                  ME = ME_list,
#  		slim = TRUE,
#                  family = poisson(),
#                  data = georgia,
#                  car_parts = cars
#                  )

## ----fig.width = 8------------------------------------------------------------
me_diag(fit, 'insurance', georgia)

## -----------------------------------------------------------------------------
mu.x <- as.matrix(fit, pars = "mu_x_true")
dim(mu.x)
head(mu.x)
mean(mu.x)

## -----------------------------------------------------------------------------
print(fit$stanfit, pars = c("mu_x_true", "car_rho_x_true", "sigma_x_true"))

## -----------------------------------------------------------------------------
x <- as.matrix(fit, pars = "x_true")
dim(x)

## -----------------------------------------------------------------------------
x.mu <- apply(x, 2, mean)
head(x.mu)

## ----fig.width = 4, fig.height = 4--------------------------------------------
rs <- resid(fit)$mean
plot(x.mu, rs)
abline(h = 0)

## ----eval = FALSE-------------------------------------------------------------
#  ME_nsp <- prep_me_data(
#    se = data.frame(insurance = georgia$insurance.se),
#    logit = TRUE
#  )
#  fit_nsp <- stan_glm(log(rate.male) ~ insurance, data = georgia, ME = ME_nsp, prior_only = TRUE)

