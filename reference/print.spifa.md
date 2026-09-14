# Print spifa Posterior Samples

Prints a fitted `spifa` object (the output of
[`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md)):
model type, formula, data dimensions, MCMC settings, and a posterior
summary table (via
[`summary.spifa`](https://ErickChacon.github.io/spifa/reference/summary.spifa.md))
– following the convention of rstan/rstanarm/brms/R2jags, which all show
actual parameter estimates by default rather than just fit metadata.

## Usage

``` r
# S3 method for class 'spifa'
print(x, ...)
```

## Arguments

- x:

  An object of class `spifa`, as returned by
  [`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md).

- ...:

  Further arguments passed to methods (currently unused).

## Value

Invisibly returns `x`.

## Author

Erick A. Chacón-Montalván

## Examples

``` r
data(ipixuna)
samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 20)
samples
#> Item factor analysis model: eifa
#> Formula: items ~ 1
#> Dimensions: 100 respondents, 10 items, 3 latent factors, 0 spatial processes
#> MCMC: 1 chain, iter = 20, thin = 1, samples = 20
#> 
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Item model parameters:
#>           mean median    sd    q10    q90 ess_bulk rhat
#> c[1]    -0.402 -0.461 0.203 -0.624 -0.140      2.8 1.49
#> c[2]    -0.499 -0.476 0.170 -0.706 -0.329      9.3 1.03
#> c[3]    -0.521 -0.569 0.200 -0.715 -0.313     10.4 0.95
#> c[4]    -0.320 -0.299 0.219 -0.649 -0.033      6.6 1.17
#> c[5]    -0.126 -0.141 0.149 -0.347  0.073     11.6 0.98
#> c[6]     0.025  0.022 0.092 -0.066  0.105     26.0 0.98
#> c[7]     0.410  0.389 0.255  0.110  0.781      7.3 1.14
#> c[8]     0.252  0.234 0.181  0.110  0.446      8.8 0.96
#> c[9]     0.613  0.643 0.161  0.406  0.823      7.1 1.16
#> c[10]    0.545  0.514 0.204  0.300  0.781      2.5 1.66
#> A[1,1]   0.835  0.909 0.194  0.508  1.061      4.1 1.22
#> A[2,1]   0.595  0.683 0.338  0.204  1.014      2.1 2.12
#> A[3,1]   0.383  0.191 0.430  0.010  0.995      2.1 2.04
#> A[4,1]   0.011 -0.200 0.379 -0.304  0.585      2.4 1.79
#> A[5,1]   0.189  0.199 0.198 -0.058  0.414      4.6 1.20
#> A[6,1]   0.256  0.188 0.239 -0.034  0.624      3.2 1.35
#> A[7,1]   0.340  0.270 0.195  0.123  0.586      3.1 1.47
#> A[8,1]  -0.249 -0.222 0.260 -0.498  0.016      2.9 1.43
#> A[9,1]   0.519  0.502 0.191  0.268  0.741      7.8 1.15
#> A[10,1]  0.359  0.382 0.219  0.040  0.643      3.0 1.46
#> A[2,2]   0.490  0.468 0.205  0.175  0.724      9.5 1.04
#> A[3,2]   0.764  0.757 0.257  0.400  1.048      8.3 1.06
#> A[4,2]   0.578  0.619 0.180  0.375  0.750      3.4 1.36
#> A[5,2]   0.637  0.577 0.381  0.244  1.128      2.1 2.12
#> A[6,2]   0.550  0.487 0.272  0.225  0.896      2.2 2.04
#> A[7,2]   0.864  0.900 0.301  0.563  1.148      9.9 1.01
#> A[8,2]   1.010  1.052 0.265  0.670  1.346      7.8 1.25
#> A[9,2]   0.762  0.783 0.315  0.281  1.133      7.5 1.17
#> A[10,2]  1.288  1.369 0.541  0.566  1.842      2.1 2.12
#> A[3,3]   0.596  0.552 0.197  0.359  0.836      2.9 1.49
#> A[4,3]   0.412  0.284 0.448 -0.046  0.957      2.3 1.91
#> A[5,3]  -0.612 -0.629 0.398 -1.043 -0.117      2.1 2.12
#> A[6,3]  -0.579 -0.622 0.373 -1.016 -0.097      2.1 2.12
#> A[7,3]  -0.349 -0.349 0.198 -0.583 -0.099      3.2 1.36
#> A[8,3]   0.284  0.258 0.193  0.074  0.530      2.5 1.69
#> A[9,3]  -0.052 -0.040 0.137 -0.222  0.092     14.0 1.11
#> A[10,3] -0.142 -0.117 0.106 -0.281 -0.031     13.2 1.13
#> 
#> ess_bulk is the bulk effective sample size; rhat is the potential
#> scale reduction factor on split chains (Rhat = 1 at convergence).
#> Use summary() for the full set of statistics (incl. ess_tail).
```
