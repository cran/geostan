#' print or plot a fitted geostan model
#'
#' @description Print a summary of model results to the R console, or plot posterior distributions of model parameters.
#'
#' @param x A fitted model object of class \code{geostan_fit}.
#' @param digits number of digits to print
#' @param probs Argument passed to \code{quantile}; which quantiles to calculate and print.
#' 
#' @param pars parameters to include; a character string (or vector) of parameter names.
#' @param plotfun Argument passed to \code{rstan::plot}. Options include histograms ("hist"), MCMC traceplots ("trace"), and density plots ("dens"). Diagnostic plots are also available such as Rhat statistics ("rhat"), effective sample size ("ess"), and MCMC autocorrelation ("ac").
#' 
#' @param fill fill color for histograms and density plots.
#'
#' @param ... additional arguments to `rstan::plot` or `rstan::print.stanfit`.
#'
#' @return
#'
#' The print methods writes text to the console to summarize model results. The plot method resturns a `ggplot` (from `rstan::plot` for stanfit objects).
#' 
#' @examples
#' data(georgia)
#' georgia$income <- georgia$income/1e3
#' 
#' fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + log(income),
#'                 centerx = TRUE,
#'                 data = georgia,
#'                 family = poisson(),
#'                 chains = 2, iter = 600) # for speed only
#'
#' 
#' # print and plot results
#' print(fit)
#' plot(fit)
#' @export
#' @md
#' @method print geostan_fit
#' @rdname print_geostan_fit
print.geostan_fit <- function(x,
                              probs = c(0.025, 0.2, 0.5, 0.8, 0.975),
                              digits = 3,
                              pars = NULL, ...) {
  pars <- c(pars, "intercept")
  cat("Spatial Model Results \n")
  cat("Formula: ")
  print(x$formula)
  if (inherits(x$slx, "formula")) {
   cat("SLX: ")
   print(x$slx)
   pars <- c(pars, "gamma")
  }
  x.pars <- c("beta", "nu", "sigma")
  if (any(x.pars %in% names(x$priors))) pars <- c(pars, names(x$priors)[grep(paste0(x.pars, collapse="$|"), names(x$priors))])   
  if(inherits(x$re$formula, "formula")) {
    cat("Partial pooling (varying intercept): ")
    print(x$re$formula)
    pars <- c(pars, "alpha_tau")
  }
  cat("Likelihood: ", x$family$family, "\n")
  cat("Link: ", x$family$link, "\n")
  
  meth <- x$spatial$method
  if (meth == "SAR") meth <- paste0(meth, " (", x$sar_type, ")")  
  cat("Spatial method: ", as.character(meth), "\n")
  if (x$spatial$method == "CAR") pars <- c(pars, "car_rho", "car_scale")
  if (x$spatial$method == "SAR") pars <- c(pars, "sar_rho", "sar_scale")  
  if (x$spatial$method == "BYM2") pars <- c(pars, "rho", "spatial_scale")
  if (x$spatial$method == "BYM") pars <- c(pars, "spatial_scale", "theta_scale")
  if (x$spatial$method == "ICAR") pars <- c(pars, "spatial_scale")  

  if (!is.null(x$diagnostic$Residual_MC)) cat("Residual Moran Coefficient: ", x$diagnostic$Residual_MC, "\n")
  cat("Observations: ", x$N, "\n")  
  if (x$ME$has_me) {
      cat("Data models (ME): ")
      cat(paste(names(x$ME$se), sep = ", ")) 
      if (x$ME$spatial_me) {
          cat("\n Data model (ME prior): CAR (auto-normal)")
      } else {
          cat("\nData model (ME prior): Student's t")
      }
  }
  if (x$spatial$method == "HS") {
    cat("\nHorseshoe global shrinkage prior: ", round(x$priors$beta_ev$global_scale, 2), "\n")
  }
  cat("\n")
  print(x$stanfit, pars = pars, digits = digits, probs = probs, ...)
}


#' @export
#' @import graphics
#' @importFrom signs signs
#' @rdname print_geostan_fit
#' @method plot geostan_fit
plot.geostan_fit <- function(x,
                             pars,
                             plotfun = "hist",
                             fill = "steelblue4",
                             ...) {
  if(missing(pars)) {
      pars <- "intercept"
      if (inherits(x$slx, "formula")) pars <- c(pars, "gamma")      
      x.pars <- c("beta", "nu", "sigma", "alpha_tau")
      if (any(x.pars %in% names(x$priors))) pars <- c(pars, names(x$priors)[grep(paste0(x.pars, collapse="|"), names(x$priors))])
      if (x$spatial$method == "SAR") pars <- c(pars, "sar_rho", "sar_scale")  
      if (x$spatial$method == "CAR") pars <- c(pars, "car_rho", "car_scale")
      if (x$spatial$method == "BYM2") pars <- c(pars, "rho", "spatial_scale")
      if (x$spatial$method == "BYM") pars <- c(pars, "spatial_scale", "theta_scale")
      if (x$spatial$method == "ICAR") pars <- c(pars, "spatial_scale")
  }
  rstan::plot(x$stanfit, pars = pars, fill = fill, plotfun = plotfun, ...) +
      scale_x_continuous(labels = signs::signs) +
      scale_y_continuous(labels = signs::signs)
}




#' Extract residuals, fitted values, or the spatial trend 
#'
#' @description Extract model residuals, fitted values, or spatial trend from a fitted \code{geostan_fit} model. 
#'
#' @param object A fitted model object of class \code{geostan_fit}.
#' 
#' @param summary Logical; should the values be summarized by their mean, standard deviation, and quantiles (\code{probs = c(.025, .2, .5, .8, .975)}) for each observation? Otherwise, a matrix containing samples from the posterior distributions is returned.
#'
#' @param rates For Poisson and Binomial models, should the fitted values be returned as rates, as opposed to raw counts? Defaults to `TRUE`; see the `Details` section for more information.
#' 
#' @param detrend For auto-normal models (CAR and SAR models with Gaussian likelihood only); if `detrend = FALSE`, the residuals will equal `y - Mu`, where `Mu` does not contain the implicit spatial trend. If `detrend = TRUE` (the default), the residuals will equal `Resid = y - Mu - Trend`. For CAR models and for SAR models of type SEM, the implicit spatial trend is `Trend = rho * C %*% (Y - Mu)` (see \link[geostan]{stan_car} or \link[geostan]{stan_sar}). For SAR models of type SLM, the spatial trend is `rho * C * y`.
#' 
#' @param trend For auto-normal models (CAR and SAR models with Gaussian likelihood only); if `trend = FALSE`, the fitted values will not include the implicit spatial trend term. Rather, the fitted values will equal `Mu`: the intercept plus any other terms explicitly included in the model formula and the `re` argument ('random effects'). If `trend = TRUE` (the default), `fitted = Mu + Trend`. For CAR models and SAR models of type SEM, the spatial trend is `rho * C (y - Mu)`. For SAR models of type SLM, the spatial trend is `rho * C * y`.
#' 
#' @param ... Not used
#'
#' @return
#'
#' By default, these methods return a `data.frame`. The column named `mean` is what most users will be looking for. These contain the fitted values (for the `fitted` method), the residuals (fitted values minus observed values, for the `resid` method), or the spatial trend (for the `spatial` method). The `mean` column is the posterior mean of each value, and the column `sd` contains the posterior standard deviation for each value. The posterior distributions are also summarized by select quantiles. 
#'
#' If `summary = FALSE` then the method returns an S-by-N matrix of MCMC samples, where S is the number of MCMC samples and N is the number of observations in the data.
#'
#' For spatial GLMs, the `spatial` method returns values on the scale of the linear predictor. Hence, for a Poisson with log link function, the spatial trend is on the log scale (in the same way that the intercept and coefficients are on the log scale).
#' 
#' @details
#'
#' When \code{rates = FALSE} and the model is Poisson or Binomial, the fitted values returned by the \code{fitted} method are the expected value of the response variable. The \code{rates} argument is used to translate count outcomes to rates by dividing by the appropriate denominator. The behavior of the `rates` argument depends on the model specification. Consider a Poisson model of disease incidence, such as the following intercept-only case:
#' ```
#' fit <- stan_glm(y ~ offset(log(E)),
#'                data = data,
#'                family = poisson())
#' ```
#' If the fitted values are extracted using `rates = FALSE`, then \code{fitted(fit)} will return the expectation of \eqn{y}. If `rates = TRUE` (the default), then \code{fitted(fit)} will return the expected value of the rate \eqn{\frac{y}{E}}.
#'
#' If a binomial model is used instead of the Poisson, then using `rates = TRUE` with the \code{fitted} method will return the expectation of \eqn{\frac{y}{N}} where \eqn{y} is the number of 'successes' and \eqn{N} is the sum of the number of 'successes' and 'failures', as in:
#' ```
#' fit <- stan_glm(cbind(successes, failures) ~ 1,
#'                data = data,
#'                family = binomial())
#' ```
#' 
#' @examples
#' sar_list <- prep_sar_data2(row = 10, col = 10, quiet = TRUE)
#' W <- sar_list$W
#' x <- sim_sar(w = W, rho = 0.8)
#' y <- sim_sar(w = W, rho = 0.7, mu = 1 + 0.5 * x)
#' dat <- data.frame(x = x, y = y)
#' 
#' fit <- stan_sar(y ~ x, data = dat, sar = sar_list,
#'                     chains = 1, iter = 800)
#'
#' # Residuals
#' r = resid(fit)
#' head(r)
#' moran_plot(r$mean, W)
#'
#' # when residuals are not detrended, they exhibit spatial autocorrelation
#' r2 = resid(fit, detrend = FALSE)
#' mc(r$mean, W)
#' mc(r2$mean, W)
#' 
#' r_mat <- resid(fit, summary = FALSE)
#' r_mean <- apply(r_mat, 2, mean)
#' head(r_mean)
#' mean(r_mean)
#' 
#' # Fitted values
#' f <- fitted(fit)
#' head(f)
#'
#' # Spatial trend
#' sp_trend  <- spatial(fit)
#' head(sp_trend)
#'
#' # another way to obtain detrended residuals 
#' sp <- spatial(fit, summary = FALSE)
#' r <- resid(fit, detrend = FALSE, summary = FALSE)
#' # detrended residuals (matrix of samples):
#' r_det <- r - sp
#' # summarize (take posterior means):
#' r_det <- apply(r_det, 2, mean)
#' r <- apply(r, 2, mean)
#' # Moran coefficients:
#' mc(r_det, W)
#' mc(r, W)
#' @export
#' @md
#' @method residuals geostan_fit
#' @rdname resid_geostan_fit
#' @export
residuals.geostan_fit <- function(object,
                                  summary = TRUE,
                                  rates = TRUE,
                                  detrend = TRUE,
                                  ...) {
    y <- object$data[,1]
    fits <- fitted(object, summary = FALSE, rates = rates, trend = detrend)
    if (rates && object$family$family == "binomial") { y <- y / (y + object$data[,2]) }
    if (rates && object$family$family == "poisson" && "offset" %in% c(colnames(object$data))) {
        log.at.risk <- object$data[, "offset"]
        at.risk <- exp( log.at.risk )
        y <- y / at.risk     
    }
    R = sweep(fits, MARGIN = 2, STATS = as.array(y), FUN = .resid)
    colnames(R) <- gsub("fitted", "residual", colnames(R))
    if (summary) return( post_summary(R) )
    return (R)  
}


#' @noRd
.resid <- function(yhat, y) y - yhat

#' @export
#' 
#' @import stats
#' @method fitted geostan_fit
#' @rdname resid_geostan_fit
fitted.geostan_fit <- function(object,
                               summary = TRUE,
                               rates = TRUE,
                               trend = TRUE,
                               ...) {
    fits <- as.matrix(object$stanfit, pars = "fitted")
    if (object$family$family == "auto_gaussian" && trend == TRUE) {
        spatial_samples <- spatial(object, summary = FALSE)
        fits <- fits + spatial_samples
    }
    if (!rates && object$family$family == "binomial") {        
        trials <- object$data[,1] + object$data[,2]
        fits <- sweep(fits, MARGIN = 2, STATS = as.array(trials), FUN = "*")
    }
    if (rates && object$family$family == "poisson" && "offset" %in% c(colnames(object$data))) {
        log.at.risk <- object$data[, "offset"]
        at.risk <- exp( log.at.risk )
        fits <- sweep(fits, MARGIN = 2, STATS = as.array(at.risk), FUN = "/")
    }
    if (summary) return( post_summary(fits) )
    return (fits)          
}

#' Extract spatial component from a fitted geostan model
#' @export
#' @rdname resid_geostan_fit
spatial <- function(object,
                    summary = TRUE,
                    ...) {
    UseMethod("spatial", object)
}

#' @export
#' @method spatial geostan_fit
#' @rdname resid_geostan_fit
spatial.geostan_fit <- function(object,
                                summary = TRUE,
                                ...) {
    if (is.na(object$spatial$par)) stop("This model does not have a spatial trend component to extract.")
    if (object$spatial$method %in% c("CAR", "SAR")) {
        spatial.samples <- extract_autoGauss_trend(object)
    } else {
        par <- as.character(object$spatial$par)
        spatial.samples <- as.matrix(object, pars = par, ...)
    }
  if (summary) {
    return ( post_summary(spatial.samples) )
  } else {
      return( spatial.samples )
  }
}

extract_autoGauss_trend <- function(object) {
    
    if (object$spatial$method == "CAR") {
        link <- object$family$link
        scale_param <- "car_scale"
        ZMP <- object$car_parts$ZMP        
        rho_name <- "car_rho"
        lag_y <- FALSE
    }
    
    if (object$spatial$method == "SAR") {
        link <- object$family$link
        scale_param <- "sar_scale"        
        ZMP <- object$sar_parts$ZMP        
        rho_name <- "sar_rho"
        lag_y <- grepl("SDLM|SLM", object$sar_type)                
    }
    
    if (object$family$family == "auto_gaussian") {

        if (lag_y) {

            C <- object$C
            y <- object$data[,1]
            Cy <- C %*% y
            rho <- as.matrix(object, pars = rho_name)
            spatial.samples <- t(sapply(1:nrow(rho), function(i) {
                as.numeric( rho[i] * Cy )
            }))
            
        } else {
            C <- object$C
            y <- object$data[,1]
            rho <- as.matrix(object, pars = rho_name)
            fits <- as.matrix(object$stanfit, pars = "fitted")
            R = sweep(fits, MARGIN = 2, STATS = as.array(y), FUN = .resid)
            spatial.samples <- t(sapply(1:nrow(rho), function(i) {
                as.numeric( rho[i] * C %*% R[i,] )
            }))
        }

        # auto-gaussan return:
        return (spatial.samples)   
    }

    # hierarchical models:
    if (ZMP == 1) {

        scale <- as.matrix(object, pars = scale_param)
        log_lambda <- as.matrix(object, pars = "log_lambda")
        spatial.samples <- sweep(log_lambda, MARGIN = 1, STATS = scale, FUN = "*")        

    } else {

        log_lambda_mu <- as.matrix(object, pars = "log_lambda_mu")
        # 'fitted, rates=TRUE' sweeps out the offset from fitted values
        lambda <- fitted(object, summary = FALSE, rates = TRUE)        
        if (link == "log") log_lambda <- log(lambda)
        if (link == "logit") log_lambda <- logit(lambda)        
        spatial.samples <- log_lambda - log_lambda_mu
        # colnames: spatial[1:N]; rownames: sample[1:S]
        row.names(spatial.samples) <- paste("sample[", 1:nrow(spatial.samples), "]")
        colnames(spatial.samples) <- paste("spatial[", 1:ncol(spatial.samples), "]")      
    }
    
    return (spatial.samples)
}



#' Extract samples from a fitted model
#'
#' @description Extract samples from the joint posterior distribution of parameters.
#'
#' @param x A fitted model object of class \code{geostan_fit}.
#'
#' @param ... Further arguments passed to \code{rstan} methods for for `as.data.frame`, `as.matrix`, or `as.array`
#'
#' @return
#'
#' A matrix, data frame, or array of MCMC samples is returned.
#' 
#' @examples
#' data(georgia)
#' A <- shape2mat(georgia, "B")
#' 
#' fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)),
#'                 C = A,
#'                 data = georgia,
#'                 family = poisson(),
#'                 chains = 1, iter = 600) # for speed only
#'
#' 
#' s <- as.matrix(fit)
#' dim(s)
#'
#' a <- as.matrix(fit, pars = "intercept")
#' dim(a)
#'
#' # Or extract the stanfit object
#' S <- fit$stanfit
#' print(S, pars = "intercept")
#' samples <- as.matrix(S)
#' dim(samples)
#' @md
#' @export
#' @rdname samples_geostan_fit
#' @method as.matrix geostan_fit
as.matrix.geostan_fit <- function(x, ...) {
  as.matrix(x$stanfit, ...)
}

#' @export
#' @rdname samples_geostan_fit
#' @method as.data.frame geostan_fit
as.data.frame.geostan_fit <- function(x, ...) {
    as.data.frame(x$stanfit, ...)
}

#' @export
#' @rdname samples_geostan_fit
#' @method as.array geostan_fit
as.array.geostan_fit <- function(x, ...){
    as.array(x$stanfit, ...)
}


#' Predict method for `geostan_fit` models
#'
#' @description Obtain predicted values from a fitted model by providing new covariate values.
#' 
#' @param object A fitted model object of class \code{geostan_fit}.
#' 
#' @param newdata A data frame in which to look for variables with which to predict. Note that if the model formula includes an offset term, `newdata` must contain the offset column (see examples below). If covariates in the model were centered using the `centerx` argument, the `predict.geostan_fit` method will automatically center the predictors in `newdata` using the values stored in `object$x_center`. If `newdata` is missing, the fitted values of the model will be returned.
#' 
#' @param alpha An N-by-1 matrix of MCMC samples for the intercept; this is provided by default. If used, note that the intercept needs to be provided on the scale of the linear predictor. This argument might be used if there is a need to incorporate the spatial trend term (as a spatially-varying intercept).
#'
#' @param center Optional vector of numeric values or a logical scalar to pass to \code{\link[base]{scale}}. Defaults to using `object$x_center`. If the model was fit using `centerx = TRUE`, then covariates were centered and their mean values are stored in `object$x_center` and the `predict` method will use them automatically to center `newdata`; if the model was fit with `centerx = FALSE`, then `object$x_center = FALSE` and `newdata` will not be centered.
#'
#' @param summary If `FALSE`, a matrix containing samples from the posterior distribution at each observation is returned. The default, `TRUE`, will summarize results by providing an estimate (mean) and credible interval (formed by taking quantiles of the MCMC samples).
#' 
#' @param type By default, results from `predict` are on the scale of the linear predictor (`type = "link")`). The alternative (`type = "response"`) is on the scale of the response variable. For example, the default return values for a Poisson model on the log scale, and using `type = "response"` will return the original scale of the outcome variable (by exponentiating the log values).
#'
#' @param add_slx Logical. If `add_slx = TRUE`, any spatially-lagged covariates that were specified through the 'slx' argument (of the model fitting function, e.g., `stan_glm`) will be added to the linear predictor. The spatial lag terms will be calculated internally using `object$C`, the spatial weights matrix used to fit the model. Hence, `newdata` must have `N = object$N` rows. Predictions from spatial lag models (SAR models of type 'SLM' and 'SDLM') always include the SLX terms (i.e., any value passed to `add_slx` will be overwritten with `TRUE`).
#' 
#' @param approx For SAR models of type 'SLM' or 'SDLM' only; use an approximation for matrix inversion? See details below.
#'
#' @param K Number of matrix powers to use with \code{approx}.
#' 
#' @param ... Not used
#'
#' @details
#'
#' The primary purpose of the predict method is to explore marginal effects of covariates. The uncertainty present in these predictions refers to uncertainty in the expected value of the model. The expectation does not include the error term of the model (nb: one expects actual observations to form a cloud of points around the expected value). By contrast, \link[geostan]{posterior_predict} returns the complete (posterior) predictive distribution of the model (the expectation plus noise).
#'
#' The model formula will be taken from `object$formula`, and then a model matrix will be created by passing `newdata` to the \link[stats]{model.frame} function (as in: \code{model.frame(object$formula, newdata)}. Parameters are taken from \code{as.matrix(object, pars = c("intercept", "beta"))}. 
#'
#' The examples illustrate how to use the function in most cases.
#'
#' Special considerations apply to models with spatially-lagged covariates and a spatially-lagged dependent variable (i.e., the 'SLM' and 'SDLM' models fit by \link[geostan]{stan_sar}). 
#' 
#' ## Spatial lag of X
#' 
#' Spatially-lagged covariates which were included via the `slx` argument will, by default, not be included in the predicted values. (The user can have greater control by manually adding the spatially-lagged covariate to the main model formula.) The `slx` term will be be included in predictions if `add_slx = TRUE` or if the fitted model is a SAR model of type 'SLM' or 'SDLM'. In either of those cases, `newdata` must have the same number of rows as were used to fit the original data.
#'
#' ## Spatial lag of Y
#' 
#' The typical 'marginal effect' interpretation of the regression coefficients does not hold for the SAR models of type 'SLM' or 'SDLM'. For details on these 'spillover effects', see LeSage and Pace (2009), LeSage (2014), and \link[geostan]{impacts}. 
#'
#' Predictions for the spatial lag model (SAR models of type 'SLM') are equal to:
#' \deqn{
#'  (I - \rho W)^{-1} X \beta
#' }
#' where \eqn{X \beta} contains the intercept and covariates. Predictions for the spatial Durbin lag model (SAR models of type 'SDLM') are equal to:
#' \deqn{
#'  (I - \rho W)^{-1} (X \beta + WX \gamma)
#' }
#' where \eqn{WX \gamma} are spatially lagged covariates multiplied by their coefficients. For this reason, the `predict.geostan_fit` method requires that `newdata` have as many rows as the original data (so that `nrow(newdata) == nrow(object$C)`); the spatial weights matrix will be taken from `object$C`. 
#'
#' The inverse of the matrix \eqn{(I - \rho W)} can be time consuming to compute (especially when iterating over MCMC samples). You can use `approx = TRUE` to approximate the inverse using a series of matrix powers. The argument \eqn{K} controls how many powers to use for the approximation. As a rule, higher values of \eqn{\rho} require larger \eqn{K} to obtain accuracy. Notice that \eqn{\rho^K} should be close to zero for the approximation to hold. For example, for \eqn{\rho = .5} a value of \eqn{K=8} may suffice (\eqn{0.5^8 = 0.004}), but larger values of \eqn{\rho} require higher values of \eqn{K}. 
#'
#' ## Generalized linear models
#' 
#' In generalized linear models (such as Poisson and Binomial models) marginal effects plots on the response scale may be sensitive to the level of other covariates in the model and to geographic location (given a spatially-varying mean value). If the model includes a spatial autocorrelation component (for example, you used a spatial CAR, SAR, or ESF model, or used the `re` argument for random effects), by default these terms will be fixed at zero for the purposes of calculating marginal effects. If you want to change this, you can introduce a varying intercept manually using the `alpha` argument.
#' 
#' @return
#'
#' If `summary = FALSE`, a matrix of samples is returned. If `summary = TRUE` (the default), a data frame is returned.
#' 
#' @examples
#' data(georgia)
#' georgia$income <- georgia$income / 1e3
#' 
#' fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + log(income),
#'                data = georgia,
#'                re = ~ GEOID,
#'                centerx = TRUE,
#'                family = poisson(),
#'                chains = 2, iter = 600) # for speed only
#'
#' # note: pop.at.risk.male=1 leads to offset of log(pop.at.risk.male)=0
#' # so that the predicted values are rates
#' newdata <- data.frame(
#'              income = seq(min(georgia$income),
#'                           max(georgia$income),
#'                            length.out = 200),
#'              pop.at.risk.male = 1)
#'
#' preds <- predict(fit, newdata, type = "response")
#' head(preds)
#' plot(preds$income,
#'      preds$mean * 10e3,
#'      type = "l",
#'      ylab = "Deaths per 10,000",
#'      xlab = "Income ($1,000s)")
#'
#' # here the predictions are rates per 10,000
#' newdata$pop.at.risk.male <- 10e3
#' preds <- predict(fit, newdata, type = "response")
#' head(preds)
#'
#' # plot range
#' y_lim <- c(min(preds$`2.5%`), max(preds$`97.5%`))
#'
#' # plot line
#' plot(preds$income,
#'     preds$mean,
#'     type = "l",    
#'     ylab = "Deaths per 10,000",
#'     xlab = "Income ($1,000s)",
#'     ylim = y_lim,
#'     axes = FALSE)
#'
#' # add shaded cred. interval
#' x <- c(preds$income, rev(preds$income))
#' y <- c(preds$`2.5%`, rev(preds$`97.5%`))
#' polygon(x = x, y = y,
#'        col = rgb(0.1, 0.2, 0.3, 0.3),
#'        border = NA)
#'
#' # add axes
#' yat = seq(0, 300, by = 20)
#' axis(2, at = yat)
#'
#' xat = seq(0, 200, by = 10)
#' axis(1, at = xat)
#'
#' # show county incomes
#' rug(georgia$income)
#' @source
#'
#' Goulard, Michael, Thibault Laurent, and Christine Thomas-Agnan (2017). About predictions in spatial autoregressive models: optimal and almost optimal strategies. *Spatial Economic Analysis* 12 (2-3): 304-325.
#' 
#' LeSage, James (2014). What Regional Scientists Need to Know about Spatial Econometrics. *The Review of Regional Science* 44: 13-32 (2014 Southern Regional Science Association Fellows Address).
#' 
#' LeSage, James, & Robert kelley Pace (2009). *Introduction to Spatial Econometrics*. Chapman and Hall/CRC.
#'
#' @importFrom Matrix Diagonal solve
#' @export
#' @md
#' 
predict.geostan_fit <- function(object,
                                newdata,
                                alpha = as.matrix(object, pars = "intercept"),                                
                                center = object$x_center,
                                summary = TRUE,
                                type = c("link", "response"),
                                add_slx = FALSE,
                                approx = FALSE,
                                K = 15,  
                                ...) {
    type <- match.arg(type)
    if (missing(newdata)) return (fitted(object, summary = summary, ...))
    f <- object$formula[-2]
    X <- as.matrix(model.matrix(f, newdata)[,-1])    
    X <- scale(X, center = center, scale = FALSE)
    O <- model.offset( model.frame(f, newdata) )
    O <- ifelse(is.null(O), 0, O)
    Beta <- as.matrix(object, pars = "beta") 
    M <- nrow(Beta)
    N <- nrow(X)
    P <- matrix(NA, nrow = M, ncol = N)
    stopifnot(nrow(alpha) == nrow(Beta))
    
    for (m in 1:M) {
        P[m,] <- O + alpha[m,] + X %*% Beta[m,]
    }

    lag_y <- FALSE
    if (object$spatial$method == "SAR") {
        add_slx <- grepl("SDEM|SDLM", object$sar_type)
        lag_y <- grepl("SLM|SDLM", object$sar_type)
    }
    
    if (add_slx) {        
        stopifnot( nrow(newdata) == nrow(object$C) )        
        Gamma <- as.matrix(object, pars = "gamma")            
        gnames <- gsub("^w.", "", colnames(Gamma))
        xnames <- colnames(X)
        idx <- na.omit( match(xnames, gnames) )
        X_tmp <- as.matrix(X[, idx])
        if (length(idx) != ncol(Gamma)) stop("Mis-match of SLX coefficient names and X column names.")
        WX <- W %*% X_tmp        
        for (m in 1:M) P[m,] <- as.numeric( P[m,] + WX %*% as.matrix(Gamma[m,]) )        
    }
    
    if (lag_y == TRUE) {
        
        if ( nrow(object$C) != N ) {
            stop("Prediction with SLM and SDLM models requires that 'newdata' is the same size as your spatial weights matrix (N rows).")
        }
        
        W <- object$C
        I <- Matrix::Diagonal(N)
        rho <- as.matrix(object, pars = "sar_rho")[,1]

        if (approx) {
            Q <- K + 1
            powers = 0:K
            Mlist <- list()
            Mlist[[1]] <- I
            Mlist[[2]] <- W
            W_k <- W    
            for (q in 3:Q) {
                W_k <- W %*% W_k
                Mlist[[q]] <- W_k
            }        
            for (m in 1:M) {
                rho_powers <- rho[ m ]^powers
                Mpowers <- lapply(seq(Q), function(j) Mlist[[j]] * rho_powers[j])
                Multiplier <- Reduce(`+`, Mpowers)
                P[ m, ] <- as.numeric( Multiplier %*% P[ m, ] )
            }           

        } else {
            for (m in 1:M) {
                Multiplier <- Matrix::solve(I - rho[m] * W) 
                P[m,] <- as.numeric( Multiplier %*% P[m,] )                    
            }
        }
        
    }
    
    if (type == "response") {
        if (object$family$link == "log") P <- exp(P)
        if (object$family$link == "logit") P <- inv_logit(P)
    }
    
    if (summary) {
        P <- post_summary(P)
        P <- cbind(newdata, P)
    }
    
    return (P)
}

#' @export
log_lik <- function(object, ...) {
        UseMethod("log_lik", object)
    }


#' @export
#' @method log_lik geostan_fit
#' @importFrom stats dt dnorm dpois
log_lik.geostan_fit <- function(object, array = FALSE, ...) {
    
    M <- object$stanfit@sim$iter - object$stanfit@sim$warmup
    K <- object$stanfit@sim$chains
    N <- object$N
    
    mu <- as.array(object, pars = "fitted")
    idx <- object$missing$y_obs_idx
    fam <- object$family$family
    
    if (fam == "student_t") {
        log_lik <- array(dim = c(M, K, N))
        y <- object$data[,'y']        
        sigma <- as.array(object, pars = "sigma")
        nu <- as.array(object, pars = "nu")        
        for (m in 1:M) {
            for (k in 1:K) {                    
                log_lik[m,k,] <- stats::dt(
                    x = (y[idx] - mu[m,k,idx])/sigma[m,k,1],
                    df = nu[m,k,1],
                    log = TRUE) - log(sigma[m,k,1])                
            }
        }
    }
    
    if (fam == "gaussian") {
        log_lik <- array(dim = c(M, K, N))        
        y <- object$data[,'y']        
        sigma <- as.array(object, pars = "sigma")
        for (m in 1:M) {
            for (k in 1:K) {
                log_lik[m,k,] <- dnorm(
                    x = y[idx],
                    mean = mu[m,k,idx],
                    sd = sigma[m,k,1],
                    log = TRUE)
            }
        }
    }
    
    if (fam == "poisson") {
        cp <- object$missing$censor_point
        if (cp == FALSE) {            
            log_lik <- array(dim = c(M, K, N))
            y <- object$data[,'y']        
            for (m in 1:M) {
                for (k in 1:K) {                    
                    log_lik[m,k,] <- dpois(
                        x = y[idx],
                        lambda = mu[m,k,idx],
                        log = TRUE)
                }
            }                               
        } else {
            N <- object$missing$n_obs + object$missing$n_mis
            log_lik <- array(dim = c(M, K, N))
            y <- object$data[,'y']
            mis_idx <- object$missing$y_mis_idx
            for (m in 1:M) {
                for (k in 1:K) {                    
                    log_lik[m,k,idx] <- dpois(
                        x = y[idx],
                        lambda = mu[m,k,idx],
                        log = TRUE)
                    log_lik[m,k,mis_idx] <- ppois(q = cp,
                                                  lambda = mu[m,k,mis_idx],
                                                  log.p = TRUE)
                }
            }            
        }       
    }
    
    if (fam == "binomial") {
        log_lik <- array(dim = c(M, K, N))
        out = model.response(model.frame(object$formula, object$data, na.action = NULL))
        y <- out[,1]
        size <- y + out[,2]        
        for (m in 1:M) {
            for (k in 1:K) {                    
                log_lik[m,k,] <- dbinom(
                    x = y[idx],
                    size = size[idx],
                    p = mu[m,k,idx],
                    log = TRUE)
            }
        }
    }

    if (fam == "auto_gaussian" & object$spatial$method == "CAR") {
        log_lik <- array(dim = c(M, K, 1))
        y <- object$data[,'y']
        sigma <- as.array(object, pars = "car_scale")
        rho <- as.array(object, pars = "car_rho")
        parts <- object$car_parts
        for (m in 1:M) {
            for (k in 1:K) {
                log_lik[m,k,1] <- car_normal_lpdf(
                    y = y,
                    mu = mu[m,k,],
                    sigma = sigma[m,k,1],
                    rho = rho[m,k,1],
                    C = parts$C,
                    D_inv = parts$Delta_inv,
                    log_det_D_inv = parts$log_det_Delta_inv,
                    lambda = parts$lambda,
                    n = parts$n
                )
            }
        }
    }

    if (fam == "auto_gaussian" & object$spatial$method == "SAR") {
        log_lik <- array(dim = c(M, K, 1))
        y <- object$data[,'y']
        sigma <- as.array(object, pars = "sar_scale")
        rho <- as.array(object, pars = "sar_rho")
        parts <- object$sar_parts
        for (m in 1:M) {
            for (k in 1:K) {
                log_lik[m,k,1] <- sar_normal_lpdf(
                    y = y,
                    mu = mu[m,k,],
                    sigma = sigma[m,k,1],
                    rho = rho[m,k,1],
                    W = parts$W,
                    lambda = parts$eigenvalues_w,
                    n = parts$n
                )
            }
        }
    }

    if (array == FALSE) log_lik <- matrix(log_lik, nrow = K * M, ncol = N)
    
    return (log_lik)
}


