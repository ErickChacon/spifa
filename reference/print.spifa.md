# Print spifa Posterior Samples

Prints a fitted `spifa` object (the output of
[`spifa`](https://ErickChacon.github.io/spifa/reference/spifa.md)):
model type, formula, data dimensions, MCMC settings, and a posterior
summary table (via
[`summary.spifa`](https://ErickChacon.github.io/spifa/reference/summary.spifa.md)).

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
#>           mean median   sd     q10     q90 ess_bulk rhat
#> c[1]    -0.412 -0.417 0.12 -0.5646 -0.2672     13.8 1.00
#> c[2]    -0.351 -0.367 0.22 -0.5662 -0.2090      8.7 1.03
#> c[3]    -0.472 -0.393 0.17 -0.6871 -0.3292      7.2 1.75
#> c[4]    -0.207 -0.227 0.14 -0.3543 -0.0058      6.3 1.14
#> c[5]    -0.196 -0.157 0.13 -0.3790 -0.0354     15.0 0.99
#> c[6]    -0.049 -0.042 0.14 -0.2137  0.1168     13.8 1.09
#> c[7]     0.262  0.216 0.16  0.1163  0.4784      3.2 1.42
#> c[8]     0.243  0.223 0.15  0.0681  0.4398     15.0 0.96
#> c[9]     0.418  0.401 0.17  0.2129  0.7170     11.2 1.03
#> c[10]    0.298  0.318 0.20  0.0456  0.5613      2.1 2.12
#> A[1,1]   0.586  0.631 0.17  0.3654  0.7787     10.2 1.40
#> A[2,1]   0.137  0.122 0.18 -0.0692  0.3631      3.0 1.43
#> A[3,1]   0.170  0.146 0.17  0.0209  0.3968      6.7 1.26
#> A[4,1]  -0.012  0.056 0.17 -0.2524  0.1550     12.2 1.00
#> A[5,1]   0.418  0.441 0.28  0.1090  0.7660      2.1 2.12
#> A[6,1]   0.249  0.286 0.18  0.0078  0.4217      3.7 1.28
#> A[7,1]   0.543  0.669 0.40  0.0290  0.9593      2.1 2.12
#> A[8,1]   0.265  0.317 0.18 -0.0184  0.4564     13.0 0.97
#> A[9,1]   0.723  0.684 0.22  0.4761  1.0157      2.1 2.12
#> A[10,1]  0.535  0.491 0.23  0.3210  0.8123      3.6 1.31
#> A[2,2]   0.755  0.719 0.20  0.5549  0.9933      2.5 1.73
#> A[3,2]   0.780  0.808 0.28  0.5300  1.1136      3.0 1.46
#> A[4,2]   0.488  0.532 0.19  0.2440  0.6985     10.5 1.06
#> A[5,2]   0.590  0.574 0.32  0.1792  0.9659      2.3 1.85
#> A[6,2]   0.540  0.574 0.23  0.1854  0.7955      2.2 1.97
#> A[7,2]   0.671  0.744 0.27  0.3736  0.9304      3.3 1.51
#> A[8,2]   0.707  0.677 0.26  0.4148  1.0296      6.6 1.13
#> A[9,2]   0.685  0.780 0.25  0.2649  0.9297     10.7 1.22
#> A[10,2]  1.038  1.090 0.44  0.3948  1.5077      4.1 1.21
#> A[3,3]   0.796  0.787 0.16  0.6303  0.9292      2.4 1.79
#> A[4,3]   0.275  0.258 0.15  0.1224  0.4499      2.5 1.69
#> A[5,3]  -0.427 -0.376 0.29 -0.8348 -0.1017      2.4 1.69
#> A[6,3]  -0.456 -0.341 0.49 -1.2241  0.0679      2.1 2.12
#> A[7,3]  -0.101 -0.108 0.20 -0.3618  0.1212      2.3 1.85
#> A[8,3]   0.623  0.589 0.42  0.0416  1.1096      2.0 2.12
#> A[9,3]  -0.031 -0.031 0.17 -0.2308  0.1459     10.4 1.05
#> A[10,3]  0.078  0.067 0.22 -0.1763  0.4005      8.5 1.35
#> 
#> ess_bulk is the bulk effective sample size; rhat is the potential
#> scale reduction factor on split chains (Rhat = 1 at convergence).
#> Use summary() for the full set of statistics (incl. ess_tail).
```
