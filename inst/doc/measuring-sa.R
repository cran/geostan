## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.align = "center",
  fig.width = 3.5,
  fig.height = 3,
  comment = "#>"
)

## ----setup, message = FALSE, warning = FALSE----------------------------------
library(geostan)
data("georgia")

## ----fig.width = 8------------------------------------------------------------
sp_diag(georgia$college, georgia, name = "College (%)")

## -----------------------------------------------------------------------------
mean(georgia$college)

## -----------------------------------------------------------------------------
weighted.mean(georgia$college, w = georgia$population)

## -----------------------------------------------------------------------------
C <- shape2mat(georgia, style = "W")
moran_plot(georgia$college, C)

## -----------------------------------------------------------------------------
mc(georgia$college, C)

## -----------------------------------------------------------------------------
moran_plot(georgia$college, shape2mat(georgia, "B"))

## -----------------------------------------------------------------------------
Li <- lisa(georgia$college, C)
head(Li)

## -----------------------------------------------------------------------------
c(mc(georgia$college, C), mean(Li$Li))

## -----------------------------------------------------------------------------
rho <- aple(georgia$ICE, C)
n <- nrow(georgia)
ess <- n_eff(rho = rho, n = n)
c(nominal_n = n, rho = rho, ESS = ess)

## -----------------------------------------------------------------------------
C <- shape2mat(georgia)
fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)), 
                data = georgia, 
                re = ~ GEOID, 
                family = poisson(), 
                C = C,
                refresh = 0 # this line silences Stan's printing
                )

## -----------------------------------------------------------------------------
print(fit)

## ----fig.width = 8------------------------------------------------------------
sp_diag(fit, georgia)

