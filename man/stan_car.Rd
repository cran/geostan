% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stan_car.R
\name{stan_car}
\alias{stan_car}
\title{Conditional autoregressive (CAR) models}
\source{
Cressie, Noel (2015 (1993)). \emph{Statistics for Spatial Data}. Wiley Classics, Revised Edition.

Cressie, Noel and Wikle, Christopher (2011). \emph{Statistics for Spatio-Temporal Data}. Wiley.

Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure. \emph{Int. J. Env. Res. and Public Health} 18 (13): 6856. DOI: 10.3390/ijerph18136856 Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.

Donegan, Connor (2021). Spatial conditional autoregressive models in Stan. \emph{OSF Preprints}. \doi{10.31219/osf.io/3ey65}.

Haining, Robert and Li, Guangquan (2020). \emph{Modelling Spatial and Spatial-Temporal Data: A Bayesian Approach}. CRC Press.
}
\usage{
stan_car(
  formula,
  slx,
  re,
  data,
  car_parts,
  C,
  family = gaussian(),
  prior = NULL,
  ME = NULL,
  centerx = FALSE,
  prior_only = FALSE,
  censor_point,
  chains = 4,
  iter = 2000,
  refresh = 500,
  keep_all = FALSE,
  pars = NULL,
  control = NULL,
  ...
)
}
\arguments{
\item{formula}{A model formula, following the R \code{\link[stats]{formula}} syntax. Binomial models can be specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.}

\item{slx}{Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally). When setting priors for \code{beta}, remember to include priors for any SLX terms.}

\item{re}{To include a varying intercept (or "random effects") term, \code{alpha_re}, specify the grouping variable here using formula syntax, as in \code{~ ID}. Then, \code{alpha_re} is a vector of parameters added to the linear predictor of the model, and:\preformatted{alpha_re ~ N(0, alpha_tau)
alpha_tau ~ Student_t(d.f., location, scale).
}

With the CAR model, any \code{alpha_re} term should be at a \emph{different} level or scale than the observations; that is, at a different scale than the autocorrelation structure of the CAR model itself.}

\item{data}{A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.}

\item{car_parts}{A list of data for the CAR model, as returned by \code{\link[geostan]{prep_car_data}}.}

\item{C}{Optional spatial connectivity matrix which will be used to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms; it will automatically be row-standardized before calculating \code{slx} terms. See \code{\link[geostan]{shape2mat}}.}

\item{family}{The likelihood function for the outcome variable. Current options are \code{auto_gaussian()}, \code{binomial(link = "logit")}, and \code{poisson(link = "log")}; if \code{family = gaussian()} is provided, it will automatically be converted to \code{auto_gaussian()}.}

\item{prior}{A named list of parameters for prior distributions (see \code{\link[geostan]{priors}}):
\describe{

\item{intercept}{The intercept is assigned a Gaussian prior distribution (see \code{\link[geostan]{normal}}}.

\item{beta}{Regression coefficients are assigned Gaussian prior distributions. Variables must follow their order of appearance in the model \code{formula}. Note that if you also use \code{slx} terms (spatially lagged covariates), and you use custom priors for \code{beta}, then you have to provide priors for the slx terms. Since slx terms are \emph{prepended} to the design matrix, the prior for the slx term will be listed first.
}

\item{car_scale}{The scale of the CAR model, \code{car_scale}. The scale is assigned a Student's t prior model (constrained to be positive).}

\item{car_rho}{The spatial autocorrelation parameter in the CAR model, \code{rho}, is assigned a uniform prior distribution. By default, the prior will be uniform over all permissible values as determined by the eigenvalues of the connectivity matrix, \code{C}. The range of permissible values for \code{rho} is automatically printed to the console by \code{\link[geostan]{prep_car_data}}.}

\item{tau}{The scale parameter for any varying intercepts (a.k.a exchangeable random effects, or partial pooling) terms. This scale parameter, \code{tau}, is assigned a Student's t prior (constrained to be positive).}

}}

\item{ME}{To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a list of data as constructed by the \code{\link[geostan]{prep_me_data}} function.}

\item{centerx}{To center predictors on their mean values, use \code{centerx = TRUE}. If the ME argument is used, the modeled covariate (i.e., latent variable), rather than the raw observations, will be centered. When using the ME argument, this is the recommended method for centering the covariates.}

\item{prior_only}{Logical value; if \code{TRUE}, draw samples only from the prior distributions of parameters.}

\item{censor_point}{Integer value indicating the maximum censored value; this argument is for modeling censored (suppressed) outcome data, typically disease case counts or deaths.}

\item{chains}{Number of MCMC chains to use.}

\item{iter}{Number of samples per chain.}

\item{refresh}{Stan will print the progress of the sampler every \code{refresh} number of samples. Set \code{refresh=0} to silence this.}

\item{keep_all}{If \code{keep_all = TRUE} then samples for all parameters in the Stan model will be kept; this is necessary if you want to do model comparison with Bayes factors and the \code{bridgesampling} package.}

\item{pars}{Optional; specify any additional parameters you'd like stored from the Stan model.}

\item{control}{A named list of parameters to control the sampler's behavior. See \code{\link[rstan]{stan}} for details.}

\item{...}{Other arguments passed to \code{\link[rstan]{sampling}}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.}
}
\value{
An object of class class \code{geostan_fit} (a list) containing:
\describe{
\item{summary}{Summaries of the main parameters of interest; a data frame.}
\item{diagnostic}{Widely Applicable Information Criteria (WAIC) with a measure of effective number of parameters (\code{eff_pars}) and mean log pointwise predictive density (\code{lpd}), and mean residual spatial autocorrelation as measured by the Moran coefficient.}
\item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
\item{data}{a data frame containing the model data}
\item{family}{the user-provided or default \code{family} argument used to fit the model}
\item{formula}{The model formula provided by the user (not including CAR component)}
\item{slx}{The \code{slx} formula}
\item{re}{A list containing \code{re}, the varying intercepts (\code{re}) formula if provided, and
\code{Data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
\item{priors}{Prior specifications.}

\item{x_center}{If covariates are centered internally (\code{centerx = TRUE}), then \code{x_center} is a numeric vector of the values on which covariates were centered.}

\item{spatial}{A data frame with the name of the spatial component parameter (either "phi" or, for auto Gaussian models, "trend") and method ("CAR")}
\item{ME}{A list indicating if the object contains an ME model; if so, the user-provided ME list is also stored here.}
\item{C}{Spatial connectivity matrix (in sparse matrix format).}
}
}
\description{
Use the CAR model as a prior on parameters, or fit data to an auto-Gaussian CAR model.
}
\details{
CAR models are discussed in Cressie and Wikle (2011, p. 184-88), Cressie (2015, Ch. 6-7), and Haining and Li (2020, p. 249-51).

The Stan code for this implementation of the CAR model first introduced in Donegan et al. (2021, supplementary material) for models of small area survey data.

Details and results depend on the \code{family} argument, as well as on the particular CAR specification chosen (see \link[geostan]{prep_car_data}).
\subsection{Auto-Gaussian}{

When \code{family = auto_gaussian()}, the CAR model is specified as follows:\preformatted{Y ~ MVGauss(Mu, Sigma)
Sigma = (I - rho C)^-1 * M * tau^2
}

where \code{Mu} is the mean vector (with intercept, covariates, etc.), \code{C} is a spatial connectivity matrix, and \code{M} is a known diagonal matrix with diagonal entries proportional to the conditional variances. \code{C} and \code{M} are provided by \code{\link[geostan]{prep_car_data}}.

The covariance matrix of the CAR model, \code{Sigma}, contains two parameters: \code{car_rho} (rho), which controls the degree of spatial autocorrelation, and the scale parameter, \code{car_scale} (tau). The range of permissible values for \code{rho} depends on the specification of \code{C} and \code{M}; for options, see \code{\link[geostan]{prep_car_data}} and Cressie and Wikle (2011, pp. 184-188).

The auto-Gaussian model contains an implicit spatial trend (i.e., autocorrelation) component which is calculated as follows (Cressie 2015, p. 564):\preformatted{trend = rho * C * (Y - Mu).
}

This term can be extracted from a fitted auto-Gaussian model using the \code{\link[geostan]{spatial}} method.

When applied to a fitted auto-Gaussian model, the \code{\link[geostan]{residuals.geostan_fit}} method returns `de-trended' residuals by default. That is,\preformatted{residual = Y - Mu - trend.
}

To obtain "raw" residuals (\code{Y - Mu}), use \code{residuals(fit, detrend = FALSE)}.
}

\subsection{Poisson}{

For \code{family = poisson()}, the model is specified as:\preformatted{Y ~ Poisson(exp(offset + lambda))
lambda ~ MVGauss(Mu, Sigma)
Sigma = (I - rho C)^-1 * M * tau^2
}

These models are most often used to calculate small area incidence rates (mortality or disease incidence rates); the user provided offset should be, then, the natural logarithm of the denominator in the rates, e.g., log-population at risk.

For Poisson models, the \code{\link[geostan]{spatial}} method returns the parameter vector \code{phi}, which is the log-risk minus the intercept and any covariates:\preformatted{phi = lambda - Mu.
}

This is the spatial autocorrelation component. This is equivalent to specifying the model as:\preformatted{Y ~ Poisson(exp(offset + Mu + phi))
phi ~ MVGauss(0, Sigma)
Sigma = (I - rho C)^-1 * M * tau^2.
}

In the Poisson CAR model, \code{phi} contains a latent spatial trend as well as additional variation around it. If you would like to extract the latent/implicit spatial trend from \code{phi}, you can do so by calculating (following Cressie 2015, p. 564):\preformatted{trend = rho * C * phi.
}
}

\subsection{Binomial}{

For \code{family = binomial()}, the model is specified as:\preformatted{Y ~ Binomial(N, theta)
logit(theta) ~ MVGauss(Mu, Sigma)
Sigma = (I - rho C)^-1 * M * tau^2
}

where outcome data \code{Y} are counts, \code{N} is the number of trials, and \code{theta} is the 'success' rate. Note that the model formula should be structured as: \code{cbind(sucesses, failures) ~ x}, such that \code{trials = successes + failures}.

For fitted Binomial models, the \code{\link[geostan]{spatial}} method will return the parameter vector \code{phi}, equivalent to:\preformatted{phi = logit(theta) - Mu.
}
}

\subsection{Spatially lagged covariates (SLX)}{

The \code{slx} argument is a convenience function for including SLX terms. For example,\preformatted{stan_glm(y ~ x1 + x2, slx = ~ x1, \\...)
}

is a shortcut for\preformatted{stan_glm(y ~ I(W \\\%*\\\% x1) + x1 + x2, \\...)
}

where \code{W} is a row-standardized spatial weights matrix (see \code{\link[geostan]{shape2mat}}). SLX terms will always be \emph{prepended} to the design matrix, as above, which is important to know when setting prior distributions for regression coefficients.

For measurement error (ME) models, the SLX argument is the only way to include spatially lagged covariates since the SLX term needs to be re-calculated on each iteration of the MCMC algorithm.
}

\subsection{Measurement error (ME) models}{

The ME models are designed for surveys with spatial sampling designs, such as the American Community Survey (ACS) estimates. With estimates, \code{x}, and their standard errors, \code{se}, the ME models have one of the the following two specifications, depending on the user input:\preformatted{x ~ Gauss(x_true, se)
x_true ~ MVGauss(mu, Sigma)
Sigma = (I - rho C)^(-1) M * tau^2
mu ~ Gauss(0, 100)
tau ~ student_t(10, 0, 40)
rho ~ uniform(lower_bound, upper_bound)
}

where the covariance matrix, \code{Sigma}, has the conditional autoregressive specification, and \code{tau} is the scale parameter. For non-spatial ME models, the following is used instead:\preformatted{x ~ Gauss(x_true, se)
x_true ~ student_t(df, mu, sigma)
df ~ gamma(3, 0.2)
mu ~ Gauss(0, 100)
sigma ~ student_t(10, 0, 40)
}

For strongly skewed variables, such census tract poverty rates, it can be advantageous to apply a logit transformation to \code{x_true} before applying the CAR or Student t prior model. When the \code{logit} argument is used, the model becomes:\preformatted{x ~ Gauss(x_true, se)
logit(x_true) ~ MVGauss(mu, Sigma)
}

and similar for the Student t model.
}

\subsection{Censored counts}{

Vital statistics systems and disease surveillance programs typically suppress case counts when they are smaller than a specific threshold value. In such cases, the observation of a censored count is not the same as a missing value; instead, you are informed that the value is an integer somewhere between zero and the threshold value. For Poisson models (\verb{family = poisson())}), you can use the \code{censor_point} argument to encode this information into your model.

Internally, \code{geostan} will keep the index values of each censored observation, and the index value of each of the fully observed outcome values. For all observed counts, the likelihood statement will be:\preformatted{p(y_i | data, model) = Poisson(y_i | fitted_i), 
}

as usual. For each censored count, the likelihood statement will equal the cumulative Poisson distribution function for values zero through the censor point:\preformatted{p(y_j | data, model) = sum_\{m=0\}^censor_point Poisson( c_m | fitted_j),
}

For example, the US Centers for Disease Control and Prevention's CDC WONDER database censors all death counts between 0 and 9. To model CDC WONDER mortality data, you could provide \code{censor_point = 9} and then the likelihood statement for censored counts would equal the summation of the Poisson probability mass function over each integer ranging from zero through 9 (inclusive), conditional on the fitted values (i.e., all model parameters). See Donegan (2021) for additional discussion, references, and Stan code.
}
}
\examples{
\donttest{
# for automatic parallel processing
#options(mc.cores = parallel::detectCores())

# model mortality rates
data(georgia)
C <- shape2mat(georgia, style = "B")
cp <- prep_car_data(C)

fit <- stan_car(deaths.male ~ offset(log(pop.at.risk.male)),
                car_parts = cp,
                data = georgia,
                family = poisson())

rstan::stan_rhat(fit$stanfit)
rstan::stan_mcse(fit$stanfit)
print(fit)
sp_diag(fit, georgia)

# censored count outcomes
sum(is.na(georgia$deaths.female))
fit <- stan_car(deaths.female ~ offset(log(pop.at.risk.female)),
                car_parts = cp,
                data = georgia,
                family = poisson(),
                censor_point = 9)

## DCAR specification (inverse-distance based)
library(sf)
A <- shape2mat(georgia, "B")
D <- sf::st_distance(sf::st_centroid(georgia))
A <- D * A
cp <- prep_car_data(A, "DCAR", k = 1)

fit <- stan_car(deaths.male ~ offset(log(pop.at.risk.male)),
               data = georgia,
               car = cp,
               family = poisson())
print(fit)
}
}
\author{
Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
}
