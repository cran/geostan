## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(eval = FALSE, echo = TRUE)

## -----------------------------------------------------------------------------
#  autonormal_file <- "autonormal.stan"

## -----------------------------------------------------------------------------
#  library(rstan)
#  library(geostan)
#  data(georgia)
#  
#  A <- shape2mat(georgia, "B")
#  car_list <- prep_car_data(A, style = "WCAR")
#  
#  # add data
#  car_list$y <- log(georgia$income / 10e3)
#  car_list$x <- matrix(log(georgia$population / 10e3), ncol = 1)
#  car_list$k <- ncol(car_list$x)
#  
#  # compile Stan model from file
#  car_model <- stan_model(autonormal_file)
#  
#  # sample from model
#  samples <- sampling(car_model, data = car_list, iter = 1e3)

## -----------------------------------------------------------------------------
#  fit <- stan_car(log(income / 10e3) ~ log(population / 10e3),
#      data = georgia, car = car_list, iter = 1e3)

## -----------------------------------------------------------------------------
#  car_poisson_file <- "car_poisson.stan"

## -----------------------------------------------------------------------------
#  library(rstan)
#  library(geostan)
#  data(georgia)
#  
#  A <- shape2mat(georgia, "B")
#  car_list <- prep_car_data(A, style = "WCAR")
#  
#  # add data
#  car_list$y <- georgia$deaths.male
#  car_list$const_offset <- log(georgia$pop.at.risk.male)
#  car_list$x <- matrix(log(georgia$income / 10e3), ncol = 1)
#  car_list$k <- ncol(car_list$x)
#  
#  # compile Stan model from file
#  car_poisson <- stan_model(car_poisson_file)
#  
#  # sample from model
#  samples <- sampling(car_poisson,
#  	data = car_list,
#  	iter = 1e3)

## -----------------------------------------------------------------------------
#  fit <- stan_car(deaths.male ~ offset(log(pop.at.risk.male)) + log(income / 10e3),
#      data = georgia, car = car_list, family = poisson(), iter = 1e3)

## -----------------------------------------------------------------------------
#  sar_model_file <- "sar_model.stan"

## -----------------------------------------------------------------------------
#  library(geostan)
#  library(rstan)
#  data(georgia)
#  
#  W <- shape2mat(georgia, "W")
#  sar_list <- prep_sar_data(W)
#  
#  # add data
#  sar_list$y <- log(georgia$income / 10e3)
#  sar_list$x <- matrix(log(georgia$population / 10e3), ncol = 1)
#  sar_list$k <- ncol(sar_list$x)
#  
#  # compile Stan model from file
#  sar_model <- stan_model(sar_model_file)
#  
#  # sample from model
#  samples <- sampling(sar_model, data = sar_list, iter = 1e3)

## -----------------------------------------------------------------------------
#  sar_dl <- prep_sar_data(W)
#  fit <- stan_sar(log(income / 10e3) ~ log(population / 10e3), data = georgia,
#      sar_parts = sar_dl, iter = 1e3)

