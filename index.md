# spifa

## Introduction

**spifa** fits item factor analysis (IFA) models for binary responses
using full Bayesian inference (Gibbs sampling with adaptive
Metropolis-Hastings), via auxiliary variables with a probit link
function. Beyond standard exploratory and confirmatory IFA, the latent
factors can be modelled as a multivariate Gaussian process to capture
spatial dependence, so spatially referenced constructs (e.g. an ideology
or socio-economic index measured at survey locations) can be mapped and
predicted at new locations.

For item factor analysis *without* spatial structure, see the
[`mirt`](https://cran.r-project.org/package=mirt) package, which
**spifa** complements rather than replaces.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("ErickChacon/spifa")
```

## Basic usage

A minimal *spatial* item factor analysis fit on the bundled `ipixuna`
dataset (an `sf` object, so a spatial Gaussian process is added
automatically – see
[`?spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md)):

``` r

library(spifa)

data(ipixuna)
samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, niter = 1000)

samples
summary(samples, burnin = 500, select = "c")
plot(samples, select = "c", burnin = 500)
```

Printing `samples` directly
([`print.spifa()`](https://ErickChacon.github.io/spifa/reference/print.spifa.md))
gives model type, dimensions, and a grouped posterior summary table at a
glance; [`summary()`](https://rdrr.io/r/base/summary.html) computes the
full set of statistics for a specific parameter block;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) gives a quick
trace + density overview.
[`plot_trace()`](https://ErickChacon.github.io/spifa/reference/plot_trace.md)/[`plot_density()`](https://ErickChacon.github.io/spifa/reference/plot_density.md)/
[`plot_interval()`](https://ErickChacon.github.io/spifa/reference/plot_interval.md)
cover trace, density, and credible-interval views individually, for any
parameter block. See
[`vignette("spifa-ipixuna")`](https://ErickChacon.github.io/spifa/articles/spifa-ipixuna.md)
for a full worked example with predictors and a theory-driven
discrimination structure.
