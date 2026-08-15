# spifa

<!-- badges: start -->
[![R-CMD-check](https://github.com/ErickChacon/spifa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ErickChacon/spifa/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Spatial Item Factor Analysis

**spifa** fits item factor analysis (IFA) models for binary responses using
full Bayesian inference (Gibbs sampling with adaptive Metropolis-Hastings),
via auxiliary variables with a probit link function. Beyond standard
exploratory and confirmatory IFA, the latent factors can be modelled as a
multivariate Gaussian process to capture spatial dependence, so spatially
referenced constructs (e.g. an ideology or socio-economic index measured at
survey locations) can be mapped and predicted at new locations.

For item factor analysis *without* spatial structure, see the
[`mirt`](https://cran.r-project.org/package=mirt) package, which **spifa**
complements rather than replaces.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ErickChacon/spifa")
```

## Usage

```r
library(spifa)

data(ipixuna)
parameters <- attr(ipixuna_wide, "parameters")
L_a <- (parameters$discrimination != 0) * 1

# data must be an sf object (its geometry becomes the spatial coordinates)
# for a spatial Gaussian process to be added to the latent factors
ipixuna_sf <- sf::st_as_sf(ipixuna_wide)
ipixuna_sf$items <- as.matrix(dplyr::select(ipixuna_wide, `Item 1`:`Item 10`))

samples <- spifa(
  items ~ x1, data = ipixuna_sf, nfactors = 2,
  niter = 1000, thin = 1, standardize = FALSE,
  constraints = list(discrimination = L_a, mgp = diag(2), resid_sd = rep(0.4^0.5, 2)))

samples_tib <- as_tibble(samples, burnin = 500)
summary(samples, burnin = 500, select = "c")
```

See `vignette("spifa-ipixuna")` for a full worked example.
