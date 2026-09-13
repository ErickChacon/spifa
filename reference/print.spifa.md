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
#> Item model parameters:
#>            mean median   sd    q10    q90 ess_bulk rhat
#> c[1]    -0.3757 -0.372 0.24 -0.660 -0.258      7.6 1.52
#> c[2]    -0.2885 -0.318 0.22 -0.523 -0.063      8.8 1.08
#> c[3]    -0.4661 -0.503 0.20 -0.691 -0.216     12.7 0.95
#> c[4]    -0.4271 -0.423 0.16 -0.629 -0.227      3.3 1.40
#> c[5]    -0.2022 -0.197 0.11 -0.354 -0.059      4.3 1.32
#> c[6]     0.0020 -0.022 0.13 -0.121  0.174     19.9 1.03
#> c[7]     0.3365  0.357 0.12  0.196  0.470      4.1 1.20
#> c[8]     0.0351  0.042 0.14 -0.132  0.231      7.7 1.10
#> c[9]     0.5487  0.525 0.19  0.326  0.807      4.5 1.19
#> c[10]    0.3935  0.398 0.16  0.184  0.611      2.7 1.60
#> A[1,1]   0.3146  0.312 0.16  0.168  0.491      3.0 1.50
#> A[2,1]  -0.3083 -0.264 0.27 -0.640 -0.082      4.6 1.20
#> A[3,1]  -0.5377 -0.520 0.22 -0.774 -0.283      2.5 1.64
#> A[4,1]  -0.8541 -0.931 0.38 -1.251 -0.337      2.4 1.69
#> A[5,1]   0.1823  0.135 0.24 -0.117  0.493      2.1 2.12
#> A[6,1]   0.0087  0.015 0.15 -0.203  0.146     12.3 1.01
#> A[7,1]  -0.3480 -0.391 0.22 -0.565 -0.024      2.3 1.85
#> A[8,1]  -0.7607 -0.865 0.33 -1.056 -0.275      2.9 1.43
#> A[9,1]  -0.3886 -0.403 0.28 -0.737 -0.016      2.1 2.12
#> A[10,1] -0.4993 -0.505 0.33 -0.884 -0.103      2.3 1.84
#> A[1,2]   0.0000  0.000 0.00  0.000  0.000       NA   NA
#> A[2,2]   1.0205  0.834 0.48  0.566  1.592      2.1 2.12
#> A[3,2]   0.6962  0.790 0.53 -0.080  1.240      2.1 2.04
#> A[4,2]   0.2304  0.252 0.24  0.015  0.487      2.2 2.12
#> A[5,2]  -0.0045 -0.025 0.21 -0.292  0.171      2.8 1.50
#> A[6,2]   0.0811  0.039 0.15 -0.111  0.268     11.3 1.00
#> A[7,2]   0.2170  0.075 0.33 -0.093  0.690      2.1 2.12
#> A[8,2]   0.3564  0.459 0.27 -0.074  0.642      2.7 1.52
#> A[9,2]   0.3876  0.330 0.19  0.225  0.590      2.6 1.65
#> A[10,2]  0.3037  0.224 0.24  0.067  0.705      3.6 1.29
#> A[1,3]   0.0000  0.000 0.00  0.000  0.000       NA   NA
#> A[2,3]   0.0000  0.000 0.00  0.000  0.000       NA   NA
#> A[3,3]   0.5757  0.575 0.22  0.326  0.792      2.2 1.97
#> A[4,3]   0.2659  0.277 0.15  0.090  0.456      4.5 1.26
#> A[5,3]   0.5892  0.599 0.27  0.275  0.853      2.2 1.97
#> A[6,3]   0.6366  0.751 0.37  0.098  1.018      2.3 1.85
#> A[7,3]   1.0916  1.226 0.43  0.478  1.551      2.1 2.12
#> A[8,3]   0.3600  0.356 0.16  0.207  0.556      6.8 1.14
#> A[9,3]   0.8707  0.929 0.30  0.477  1.163      2.5 1.62
#> A[10,3]  0.9517  1.011 0.40  0.378  1.433      2.5 1.59
#> 
#> ess_bulk is the bulk effective sample size; rhat is the potential
#> scale reduction factor on split chains (Rhat = 1 at convergence).
#> Use summary() for the full set of statistics (incl. ess_tail).
```
