
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" align="right" width="160" />

## geostan: Bayesian spatial analysis

The [**geostan**](https://connordonegan.github.io/geostan/) R package
supports a complete spatial analysis workflow with Bayesian models for
areal data, including a suite of functions for visualizing spatial data
and model results. For demonstrations and discussion, see the package
[help
pages](https://connordonegan.github.io/geostan/reference/index.html) and
[vignettes](https://connordonegan.github.io/geostan/articles/index.html)
on spatial autocorrelation and spatial measurement error models.

The package is designed primarily to support public health research with
spatial data; see the
[**surveil**](https://connordonegan.github.io/surveil/) R package for
time series analysis of public health surveillance data.

**geostan** is an interface to [**Stan**](https://mc-stan.org), a
state-of-the-art platform for Bayesian inference.

### Disease mapping and spatial regression

Statistical models for data recorded across areal units like states,
counties, or census tracts.

### Observational uncertainty

Incorporate information on data reliability, such as standard errors of
American Community Survey estimates, into any **geostan** model.

### Censored observations

Vital statistics and disease surveillance systems like CDC Wonder censor
case counts that fall below a threshold number; **geostan** can model
disease or mortality risk with censored observations.

### Spatial analysis tools

Tools for visualizing and measuring spatial autocorrelation and map
patterns, for exploratory analysis and model diagnostics. Visual
diagnostics also support the evaluation of survey data quality and
observational error models.

### The RStan ecosystem

Compatible with a suite of high-quality R packages for Bayesian
inference and model evaluation.

### Custom spatial models

Tools for building custom spatial models in
[Stan](https://mc-stan.org/).

## Installation

Install **geostan** from CRAN using:

``` r
install.packages("geostan")
```

## Usage

Load the package and the `georgia` county mortality data set (ages
55-64, years 2014-2018):

``` r
library(geostan)
data(georgia)
```

The `sp_diag` function provides visual summaries of spatial data,
including a histogram, Moran scatter plot, and map:

``` r
A <- shape2mat(georgia, style = "B")
sp_diag(georgia$rate.female, georgia, w = A)
#> 3 NA values found in x; they will be dropped from the data before creating the Moran plot. If matrix w was row-standardized, it no longer is. You may want to use a binary connectivity matrix using style = 'B' in shape2mat.
#> Warning: Removed 3 rows containing non-finite values (stat_bin).
```

<img src="man/figures/README-unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

There are three censored observations in the `georgia` female mortality
data, which means there were 9 or fewer deaths in those counties. The
following code fits a spatial conditional autoregressive (CAR) model to
female county mortality data. By using the `censor_point` argument we
include our information on the censored observations to obtain results
for all counties:

``` r
cars <- prep_car_data(A)
#> Range of permissible rho values:  -1.661134 1
fit <- stan_car(deaths.female ~ offset(log(pop.at.risk.female)),
                censor_point = 9,
        data = georgia,
        car_parts = cars,
        family = poisson(),
        cores = 4, # for multi-core processing
        refresh = 0) # to silence some printing
#> 
#> *Setting prior parameters for intercept
#> Distribution: normal
#>   location scale
#> 1     -4.7     5
#> 
#> *Setting prior for CAR scale parameter (car_scale)
#> Distribution: student_t
#>   df location scale
#> 1 10        0     3
#> 
#> *Setting prior for CAR spatial autocorrelation parameter (sar_rho)
#> Distribution: uniform
#>   lower upper
#> 1  -1.7     1
```

Passing a fitted model to the `sp_diag` function will return a set of
diagnostics for spatial models:

``` r
sp_diag(fit, georgia, w = A)
#> 3 NA values found in x; they will be dropped from the data before creating the Moran plot. If matrix w was row-standardized, it no longer is. You may want to use a binary connectivity matrix using style = 'B' in shape2mat.
#> Warning: Removed 3 rows containing missing values (geom_pointrange).
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

The `print` method returns a summary of the probability distributions
for model parameters, as well as Markov chain Monte Carlo (MCMC)
diagnostics from Stan (Monte Carlo standard errors of the mean
`se_mean`, effective sample size `n_eff`, and the R-hat statistic
`Rhat`):

``` r
print(fit)
#> Spatial Model Results 
#> Formula: deaths.female ~ offset(log(pop.at.risk.female))
#> Spatial method (outcome):  CAR 
#> Likelihood function:  poisson 
#> Link function:  log 
#> Residual Moran Coefficient:  0.0017505 
#> WAIC:  1290.39 
#> Observations:  159 
#> Data models (ME): none
#> Inference for Stan model: foundation.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>             mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff  Rhat
#> intercept -4.673   0.002 0.089 -4.851 -4.720 -4.675 -4.630 -4.488  2465 1.001
#> car_rho    0.926   0.001 0.057  0.783  0.898  0.939  0.969  0.995  2896 1.001
#> car_scale  0.456   0.001 0.036  0.391  0.432  0.454  0.479  0.532  3110 1.001
#> 
#> Samples were drawn using NUTS(diag_e) at Mon Sep 19 15:23:26 2022.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

More demonstrations can be found in the package [help
pages](https://connordonegan.github.io/geostan/reference/index.html) and
[vignettes](https://connordonegan.github.io/geostan/articles/index.html).