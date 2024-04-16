
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" align="right" width="160" />

## geostan: Bayesian spatial analysis

The [**geostan**](https://connordonegan.github.io/geostan/) R package
supports a complete spatial analysis workflow with Bayesian models for
areal data, including a suite of functions for visualizing spatial data
and model results. **geostan** models were built using
[**Stan**](https://mc-stan.org), a state-of-the-art platform for
Bayesian modeling. The package is designed partly for public health
research with spatial data, for which it complements the
[**surveil**](https://connordonegan.github.io/surveil/) R package for
time series analysis of public health surveillance data.

Features include:

  - **Disease mapping and spatial regression** Statistical models for
    data recorded across areal units like states, counties, or census
    tracts.
  - **Spatial analysis tools** Tools for visualizing and measuring
    spatial autocorrelation and map patterns, for exploratory analysis
    and model diagnostics.  
  - **Observational uncertainty** Incorporate information on data
    reliability, such as standard errors of American Community Survey
    estimates, into any **geostan** model.
  - **Missing and Censored observations** Vital statistics and disease
    surveillance systems like CDC Wonder censor case counts that fall
    below a threshold number; **geostan** can model disease or mortality
    risk for small areas with censored observations or with missing
    observations.
  - **The RStan ecosystem** Interfaces easily with many high-quality R
    packages for Bayesian modeling.
  - **Custom spatial models** Tools for building custom spatial models
    in [Stan](https://mc-stan.org/).

## Installation

Using your R console, you can install **geostan** from CRAN:

``` r
install.packages("geostan")
```

You can install the latest version from the package github repository:

``` r
if (!require('devtools')) install.packages('devtools')
devtools::install_github("connordonegan/geostan")
```

If you are using Windows, you may need to install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) first. To
install Rtools:

1.  Visit the Rtools site:
    <https://cran.r-project.org/bin/windows/Rtools/>
2.  Select the version that corresponds to the version of R that you
    have installed (e.g., R 4.3).
3.  After selecting the correct version, look for the “Install Rtools”
    section (just below the introductory text) and click on the
    “installer” to download it. For example, for Rtools43 (for R
    version 4.3), click on “<span style="color:blue">Rtools43
    installer</span>.”
4.  Go to the `.exe` file you just downloaded and double-click to begin
    installation of Rtools.

## Support

All functions and methods are documented (with examples) on the website
[reference](https://connordonegan.github.io/geostan/reference/index.html)
page. See the package
[vignettes](https://connordonegan.github.io/geostan/articles/index.html)
for more on exploratory spatial analysis, spatial measurement error
models, spatial regression with raster layers, and building custom
spatial model in Stan.

To ask questions, report a bug, or discuss ideas for improvements or new
features please visit the
[issues](https://github.com/ConnorDonegan/geostan/issues) page, start a
[discussion](https://github.com/ConnorDonegan/geostan/discussions), or
submit a [pull request](https://github.com/ConnorDonegan/geostan/pulls).

## Usage

Load the package and the `georgia` county mortality data set:

``` r
library(geostan)
#> This is geostan version 0.6.0
data(georgia)
```

This has county population and mortality data by sex for ages 55-64, and
for the period 2014-2018. As is common for public access data, some of
the observations missing because the CDC has censored them.

The `sp_diag` function provides visual summaries of spatial data,
including a histogram, Moran scatter plot, and map. Here is a visual
summary of crude female mortality rates (as deaths per 10,000):

``` r
A <- shape2mat(georgia, style = "B")
mortality_rate <- georgia$rate.female * 10e3
sp_diag(mortality_rate, georgia, w = A)
#> 3 NA values found in x will be dropped from data x and matrix w
#> Warning: Removed 3 rows containing non-finite outside the scale range
#> (`stat_bin()`).
```

<img src="man/figures/README-unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

The following code fits a spatial conditional autoregressive (CAR) model
to female county mortality data. These models are used for estimating
disease risk in small areas like counties, and for analyzing covariation
of health outcomes with other area qualities. The R syntax for fitting
the models is similar to using `lm` or `glm`. We provide the population
at risk (the denominator for mortality rates) as an offset term, using
the log-transform. In this case, three of the observations are missing
because they have been censored; per CDC criteria, this means that there
were 9 or fewer deaths in those counties. By using the `censor_point`
argument and setting it to `censor_point = 9`, the model will account
for the censoring process when providing estimates of the mortality
rates:

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
#> 3 NA values identified in the outcome variable
#> Found in rows: 55, 126, 157
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
#> *Setting prior for CAR spatial autocorrelation parameter (car_rho)
#> Distribution: uniform
#>   lower upper
#> 1  -1.7     1
```

Passing a fitted model to the `sp_diag` function will return a set of
diagnostics for spatial models:

``` r
sp_diag(fit, georgia, w = A)
#> Using sp_diag(y, shape, rates = TRUE, ...). To examine data as (unstandardized) counts, use rates = FALSE.
#> 3 NA values found in x will be dropped from data x and matrix w
#> Warning: Removed 3 rows containing missing values or values outside the scale
#> range (`geom_pointrange()`).
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
#> Residual Moran Coefficient:  0.002382 
#> WAIC:  1228.61 
#> Observations:  156 
#> Data models (ME): none
#> Inference for Stan model: foundation.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>             mean se_mean    sd   2.5%    20%    50%    80%  97.5% n_eff  Rhat
#> intercept -4.677   0.003 0.099 -4.855 -4.734 -4.676 -4.615 -4.498  1384 1.002
#> car_rho    0.926   0.001 0.057  0.781  0.886  0.940  0.974  0.996  2921 1.001
#> car_scale  0.457   0.001 0.034  0.395  0.428  0.455  0.485  0.526  3929 1.000
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Apr 16 08:16:01 2024.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Applying the `fitted` method to the fitted model will return the fitted
values from the model - in this case, the fitted values are the
estimates of the county mortality rates. Multiplying them by 10,000
gives mortality rate per 10,000 at risk:

``` r
mortality_est <- fitted(fit) * 10e3
county_name <- georgia$NAME
head( cbind(county_name, mortality_est) )
#>           county_name      mean        sd      2.5%       20%       50%
#> fitted[1]       Crisp 101.70572  9.498562  83.97063  93.70448 101.24639
#> fitted[2]     Candler 137.39050 16.606716 106.86044 123.51343 136.32448
#> fitted[3]      Barrow  94.24360  6.351406  82.25888  88.85746  94.05780
#> fitted[4]      DeKalb  59.76307  1.544278  56.72223  58.42850  59.78620
#> fitted[5]    Columbia  53.31970  3.321858  47.13175  50.44742  53.30227
#> fitted[6]        Cobb  54.12199  1.467704  51.33636  52.85897  54.10105
#>                 80%     97.5%
#> fitted[1] 109.56469 121.40981
#> fitted[2] 150.82533 173.49393
#> fitted[3]  99.61073 107.25586
#> fitted[4]  61.08625  62.79658
#> fitted[5]  56.12651  59.92468
#> fitted[6]  55.34254  57.02782
```

The mortality estimates are stored in the column named “mean”, and the
limits of the 95% credible interval are found in the columns “2.5%” and
“97.5%”.

Details and demonstrations can be found in the package [help
pages](https://connordonegan.github.io/geostan/reference/index.html) and
[vignettes](https://connordonegan.github.io/geostan/articles/index.html).

## Citing geostan

If you use geostan in published work, please include a citation.

Donegan, Connor (2022) “geostan: An R package for Bayesian spatial
analysis” *The Journal of Open Source Software*. 7, no. 79: 4716.
<https://doi.org/10.21105/joss.04716>.

[![DOI](https://joss.theoj.org/papers/10.21105/joss.04716/status.svg)](https://doi.org/10.21105/joss.04716)

    @Article{,
      title = {{geostan}: An {R} package for {B}ayesian spatial analysis},
      author = {Connor Donegan},
      journal = {The Journal of Open Source Software},
      year = {2022},
      volume = {7},
      number = {79},
      pages = {4716},
      doi = {10.21105/joss.04716},
    }
