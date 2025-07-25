# geostan 0.8.2

This fixes #25, an issue with the `geostan::spatial` method for GLMs, and also includes some minor changes to package documentation.

# geostan 0.8.1

This is a small patch to ensure that the package can install on older versions of R. 

The patch changes one line in the 'Spatial analysis with geostan' vignette that caused the package to fail on installation for old versions of windows and mac os. 

# geostan 0.8.0

New Features include spatial econometric models and better sampling more hierarchical spatial autoregressive models:

 - geostan now provides more support for spatial econometric models. Only the most used models are included: spatial error (SEM), spatially-lagged dependent variable (SLM), and the 'Durbinized' versions of each (SDEM, SDLM). Methods for calculating impacts are and predicted values for SLM/SDLM are included. Users can also draw from the posterior predictive distributions of these models: see `geostan::stan_sar` and `geostan::impacts`.
 - `geostan::sim_sar` can now simulate draws from the spatial lag model as well as the spatial error model.
 - The CAR and SAR models, when placed in a Poisson or binomial model for count data, can sample poorly when the data is sparse. This can usually be fixed using the new `zmp` option which switches CAR/SAR models to a zero-mean parameterization. Details are in the vignette on building custom spatial models.
 
Bug fix:

 - A bug was introduced to v0.7.0 which added a small bias to coefficient estimates which only appeared when the `slx` option was used. This has been fixed.

Other updates:

 - A set of Monte Carlo analyses have been added to the package tests, verifying that the estimates match expectations. This is also part of a process of moving away from the 'testthat' package. 
 - The package documentation and vignettes have also been improved a little bit.


# geostan 0.7.0

V0.7.0 includes various adjustments to speed things up, and the DIC is now provided for model comparison (in addition to WAIC).

New features:

 - Deviance information criteria (DIC). WAIC is considered to be more robust than DIC, but DIC remains common in spatial statistics because (unlike WAIC) it does not assume independent observations.

Changes in the background that should improve the user experience:

 - Samples for the pointwise log likelihood are no longer collected in the Stan models. Instead, they are calculated in R when needed for calculating waic or dic. This should improve sampling speed and lower memory and storage use.
 - The QR decomposition is now implemented (automatically) for models with covariates. (Note that this cannot be used for measurement error models.) The benefits of the QR decomposition in terms of sampling speed and efficiency are greatest when covariates are centered (which you can do using the `centerx` argument). 
 - A change recommended by Roger Bivand takes advantage of an improvement in `spdep`'s creation of neighbors objects (https://github.com/ConnorDonegan/geostan/issues/19). This will speed up the `shape2mat` function in some cases.

Other changes:

 - Stan code for the CAR and SAR models has been simplified a bit. The changes are reflected in the vignette on custom spatial models. The output of `prep_car_data` and `prep_sar_data` have changed somewhat, but the user workflow is the same.

# geostan 0.6.2

geostan was removed from CRAN for a moment due to an issue with the StanHeaders R package. This should be resolved now. This release puts geostan back on CRAN with only minimal internal changes to geostan. 

# geostan 0.6.1

Two updates:

1. Multiple changes related to spatial connectivity matricies, including a new vignette.
2. A change to `geostan::predict`.

There are three updates related to spatial connectivity matrices:

 - There is a new vignette on spatial connectivity matrices (see `browseVignettes('geostan')`), written for new users.
 - Visualizing spatial neighbors: `geostan::edges` can now return a simple features object; this can be used to visualize (map) the graph structure of the spatial connectivity matrix. There is an example in the new vignette.
 - Changes to `geostan::shape2mat`: an option for k-nearest neighbors has been added, the `queen` argument is being replaced by `method`, and the function now prints a summary of the matrix to the console (using the new `geostasn::n_nbs` function)

There was one change to the `geostan::predict` method:

  - this method is for getting predicted values from a fitted model, given user-provided covariate values. Previously, the function used a point estimate of the intercept when calculating the predicted values. Now, the method uses samples from the posterior distribution of the intercept, just like it does for other parameters. The user (still) has the option to provide a matrix of samples for the intercept term; this is in case one wishes to incorporate the spatial trend term into the intercept (as a spatially-varying intercept).

# geostan 0.6.0

Updates:

1. Missing outcome data is now allowed in most models
2. A bug in `prep_icar_data` has been fixed

The model fitting functions (`stan_glm`, `stan_car`, etc.) now allow for missing data in the outcome variable. This is explained in the `geostan::stan_glm` documentation, next to the discussion of handling censored observations. When missing observations are present, there will (only) be a warning issued. This functionality is available for any GLM (`stan_glm`), any ESF model (`stan_esf`), and any model for count data (Poisson and binomial models including CAR and SAR models). The only models for which this functionality is not currently available are CAR and SAR models that are being been fit to continuous outcome variables.

The `prep_icar_data` function, which is used inside `stan_icar`, did not have the expected behavior in all cases - this has been fixed thanks to this [pull request](https://github.com/ConnorDonegan/geostan/pull/18).

The package home page now has instructions for installing from github using `devtools::install_github` https://connordonegan.github.io/geostan/

# geostan 0.5.4

Minor updates to the vignettees and documentation, also re-compiled geostan models using the latest StanHeaders (fixing an error on CRAN).

# geostan 0.5.3

## Minor changes

The `gamma` function (which is available to help set prior distributions) has been renamed to `geostan::gamma2` to avoid conflict with `base::gamma`. 

Some code for `geostan::stan_car` was cleaned up to avoid sending duplicate variables to the Stan model when a spatial ME (measurement error) model was used: https://github.com/ConnorDonegan/geostan/issues/17. This should not change any functionality and there is no reason to suspect that results were ever impacted by the duplicate variables. 

# geostan 0.5.2

This release was built using rstan 2.26.23, which incorporates Stan's new syntax for declaring arrays. Some models seems to run a little bit faster, but otherwise there are no changes that users should notice.

The warnings issued about the sp package can be ignored; these are due to geostan's dependence on spdep, which imports sp but does not use any of the deprecated functions. 

A new vignette shows how to implement some of geostan's spatial models directly in Stan, using the custom Stan functions that make the CAR and SAR models sample quickly, and using some geostan functions that make the data cleaning part easy. 

# geostan 0.5.1

### Minor fixes

This release fixes some issues that were introduced with the `slim` and `drop` arguments (in v0.5.0).

# geostan 0.5.0

## New additions

The package now provides some support for spatial regression with raster data, including for layers with hundreds of thousands of observations (possibly more, depending on one's computational resources). Two new additions make this possible.

 1. `slim = TRUE` The model fitting functions (`stan_glm`, `stan_car`, `stan_sar`, `stan_esf`, `stan_icar`) now provide the option to trim down the parameters for which MCMC samples are collected. For large N and/or many N-length vectors of parameters, this option can speed up sampling considerably and reduce memory usage. The new `drop` argument provides users control over which parameter vectors will be ignored. This functionality may be helpful for any number of purposes, including modeling large data sets, measurement error models, and Monte Carlo studies.
 2. `prep_sar_data2` and `prep_car_data2` These two functions can quickly prepare required data for SAR and CAR models when using raster layers (observations on a regularly spaced grid). The standard and more generally applicable functions `prep_car_data` and `prep_sar_data` are limited in terms of the size of spatial weights matrices they can handle.

These new functions are discussed in a new vignette titled "Raster regression."

## Minor changes

The PDF documentation has been improved---previously, multi-line equations were not rendered properly. Now they render correctly, and a mistake in the description of Binomial CAR models has been corrected.

# geostan 0.4.1

## Minor changes

 - The recommended citation for the software has been updated since the software has gone through peer-review in *The Journal of Open Source Software*. Many thanks to the two peer reviewers of the project, Chris Jochem and Virgilio Gómez Rubio. The following changes were introduced following Chris J.'s recommendations.
 - The spatial diagonstic function (`sp_diag`) will now take a spatial connectivity matrix from the fitted model object provided by the user. This way the matrix will be the same one that was used to fit the model. (All of the model fitting functions have been updated to support this functionality.) 
 - The documentation of the methods for fitted models (`residuals`, `fitted`, `spatial`, etc.) were previously packed into one page. Now, the documentation is spread over a few pages and the methods are grouped together in a more reasonable fashion. 

# geostan 0.4.0

## New Additions

### SAR models

The simultaneously-specified spatial autoregressive (SAR) model---referred to as the spatial error model (SEM) in the spatial econometrics literature---has been implemented. The SAR model can be applied directly to continuous data (as the likelihood function) or it can be used as prior model for spatially autocorrelated parameters. Details are provided on the documentation page for the `stan_sar` function.

## Minor changes

 - Previously, when getting fitted values from an auto-normal model (i.e., the CAR model with `family = auto_gaussian()`) the fitted values did not include the implicit spatial trend. Now, the `fitted.geostan_fit` method will return the fitted values with the implicit spatial trend; this is consistent with the behavior of `residuals.geostan_fit`, which has an option to `detrend` the residuals. This applies to the SAR and CAR auto-normal specifications. For details, see the documentation pages for `stan_car` and `stan_sar`.

 - The documentation for the models (`stan_glm`, `stan_car`, `stan_esf`, `stan_icar`, `stan_sar`) now uses Latex to typeset the model equations.

# geostan 0.3.0

## New additions

 - New exploratory spatial data analysis functions have been added: the Geary Ratio (GR) and the local Geary's C. These complement the Moran coefficient and local Moran's I. 
 - The vignette on spatial autocorrelation has been updated and expanded, including with a short discussion of exploratory spatial data analysis (ESDA). 
 - The vignette on spatial measurement error models/working with ACS data has been completely re-written.

## Minor changes

 - geostan models can now be used with the bridgesampling package for model comparison with Bayes factors (e.g., use `bridge_sampler(geostan_fit$stanfit)`). By default, geostan only collects MCMC samples for parameters that are expected to be of some interest for users. To become compatible with bridgesampling, the `keep_all` argument was added to all of the model fitting functions. For important background and details see the bridgesampling package documentation and vignettes on [CRAN](https://CRAN.R-project.org/package=bridgesampling).
 - stan_car now has an option to provide the connectivity matrix C, which is used to calculate spatial-lag of X (SLX) terms and residual spatial autocorrelation. Previously, there was no option to provide this matrix, as it was taken from the car_parts argument. However, that choice is only appropriate when the WCAR specification is used. Now, if C is missing and the WCAR specification has not been used a warning will appear.
 - Previously, the `lisa` function would automatically center and scale the variate before computing local Moran's I. Now, the variate will be centered and scaled by default but the user has the option to turn the scaling off (so the variate will be centered, but not divided by its standard deviation). This function also row-standardized the spatial weights matrix automatically, but there was no reason why. That's not done anymore.

# geostan 0.2.1

## Minor changes

The distance-based CAR models that are prepared by the `prep_car_data` function have changed slightly. The conditional variances were previously a function of the sum of neighboring inverse distances (in keeping with the specification of the connectivity matrix); this can lead to very skewed frequency distributions of the conditional variances. Now, the conditional variances are equal to the inverse of the number of neighboring sites. This is in keeping with the more common CAR model specifications.

# geostan 0.2.0

## Major changes

### Models for censored disease and mortality data

geostan now supports Poisson models with censored count data, a common problem in public health research where small area disease and mortality counts are censored below a threshold value. Model for censored outcome data can now be implemented using the `censor_point` argument found in all of the model fitting functions (stan_glm, stan_car, stan_esf, stan_icar).

### Measurement error models improved

The measurement error models have been updated in three important respects:

  - There is now a prep_me_data function which must be used to create the list of data for the ME models. See `?prep_me_data`.
  - For covariates that are proportions or rates, the ME models now have an option for using a logit transformation on the variable. Again, see `?prep_me_data` for usage.
  - Previously, when using `stan_car`, ME models automatically employed the CAR model as a prior for the modeled covariates. That has changed, so that the default behavior for the ME models is the same across all `stan_*` models (CAR, GLM, ESF, ICAR). 

The second change listed above is particularly useful for variables that are highly skewed, such as the poverty rate. To determine whether a transformation should be considered, it can be helpful to evaluate results of the ME model (with the untransformed covariate) using the `me_diag` function. The logit transform is done on the 'latent' (modeled) variable, not the raw covariate. This transformation cannot be applied to the raw data by the user because that would require the standard errors of covariate estimates (e.g., ACS standard errors) to be adjusted for the transformation.

## Minor changes

### A predict method for marginal effects

A `predict` method has been introduced for fitted geostan models; this is designed for calculating marginal effects. Fitted values of the model are still returned using `fitted` and the posterior predictive distribution is still accessible via `posterior_predict`.

### Centering covariates with measurement error models

The `centerx` argument has been updated to handle measurement error models for covariates. The centering now happens inside the Stan model so that the means of the modeled covariates (latent variables) are used instead of the raw data mean. 

# geostan 0.1.1

## Minor changes

  - The stan files for the CAR model have been combined with the 'foundation.stan' file, which compresses the file size considerably.
  - The vignette on spatial autocorrelation has also been updated to include model diagnostics.
  - A new example has been added to the stan_car documentation.

# geostan 0.1.0

geostan's first release.

