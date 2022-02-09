## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      eval = TRUE, 
                      fig.align = "center",
                      fig.width = 3.5,
                      fig.height = 3
                      )

## ----message = FALSE, warning = FALSE-----------------------------------------
library(geostan)
library(ggplot2)
theme_set(theme_classic())
data(georgia)

## ----fig.width = 7------------------------------------------------------------
sp_diag(georgia$ICE.se, georgia, name = "SE(ICE)")

## -----------------------------------------------------------------------------
c(sd.ice <- sd(georgia$ICE))
c(mad.ice <- mad(georgia$ICE))
scaled_se <- georgia$ICE.se / mad.ice
ggplot() +
  geom_histogram(aes(scaled_se),
                 col = 'gray50',
                binwidth = 0.05
                 )

## ----eval = FALSE-------------------------------------------------------------
#  fit <- stan_glm(log(rate.male) ~ ICE, data = georgia)

## -----------------------------------------------------------------------------
# use binary weights matrix for prep_car_data
C <- shape2mat(georgia, style = "B")
cp <- prep_car_data(C, style = "WCAR")
ME <- prep_me_data(
  se = data.frame(ICE = georgia$ICE.se),
  car_parts = cp,
  bounds = c(-1, 1)
)

## -----------------------------------------------------------------------------
ME <- prep_me_data(
  se = data.frame(ICE = georgia$ICE.se),
  car_parts = cp,
  bounds = c(-1, 1),
  prior = list(location = normal(0, 0.5),
  	       scale = student_t(df = 10, location = 0, scale = 1))
)

## -----------------------------------------------------------------------------
fit <- stan_glm(log(rate.male) ~ ICE, data = georgia, ME = ME, prior_only = TRUE)

## ----fig.width = 7.5----------------------------------------------------------
# (math symbols in figure labels may not be visible on html versions of this vignette)
me_diag(fit, 'ICE', georgia)

## ----fig.width = 7.5----------------------------------------------------------
me_diag(fit, 'ICE', georgia, index = 5)

## -----------------------------------------------------------------------------
delta <- me_diag(fit, 'ICE', georgia, plot = FALSE)$delta_data
head(delta)

## -----------------------------------------------------------------------------
georgia[c(91, 105, 90, 39), c("NAME", "population", "white", "black", "hisp", "ai", "ICE", "ICE.se", "college", "college.se")]

## -----------------------------------------------------------------------------
mu.x <- as.matrix(fit, pars = "mu_x_true")
dim(mu.x)
mean(mu.x)

## -----------------------------------------------------------------------------
print(fit$stanfit, pars = c("mu_x_true", "car_rho_x_true", "sigma_x_true"))

## -----------------------------------------------------------------------------
x <- as.matrix(fit, pars = "x_true")
dim(x)

## -----------------------------------------------------------------------------
x.mu <- apply(x, 2, mean)
head(x.mu)

## ----eval = FALSE-------------------------------------------------------------
#  ME_nsp <- prep_me_data(
#    se = data.frame(ICE = georgia$ICE.se),
#    bounds = c(-1, 1)
#  )
#  fit_nsp <- stan_glm(log(rate.male) ~ ICE, data = georgia, ME = ME_nsp, prior_only = TRUE)

## ----eval = TRUE--------------------------------------------------------------
fit_2 <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + ICE, 
                  re = ~ NAME,   
                  data = georgia, 
                  ME = ME, 
                  family = poisson(), 
                  refresh = 0)

## -----------------------------------------------------------------------------
print(fit_2)

## ----fig.width = 4.5----------------------------------------------------------
ggplot(georgia) +
  geom_point(aes(ICE, log(rate.male), col = ICE.se),
             shape = 6,
             lwd = 2) +
  labs(x = "ICE Estimate", y = "Crude log mortality") +
  scale_colour_gradient(low = "white", high = "darkred", name = "SE(ICE)") +
  theme(panel.background = element_rect(fill = 'gray20'),
        plot.background = element_rect(fill = 'gray80'),
        legend.background = element_rect(fill = 'gray80')
        )

## ----fig.width = 4.5----------------------------------------------------------
eta <- fitted(fit_2)
x <- as.matrix(fit_2, pars = "x_true")
x.mean <- apply(x, 2, mean)
x.sd <- apply(x, 2, sd)

ggplot(georgia) +
  geom_point(aes(x.mean, log(eta$mean), col = x.sd),
             shape = 6,
             lwd = 2) +
  labs(x = "ICE", y = "Log mortality") +
  scale_colour_gradient(low = "white", high = "darkred", name = "SD(ICE)") +
  theme(panel.background = element_rect(fill = 'gray20'),
        plot.background = element_rect(fill = 'gray80'),
        legend.background = element_rect(fill = 'gray80')
        )

