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
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Item model parameters:
#>             mean median    sd    q10    q90 ess_bulk rhat
#> c[1]    -0.61827 -0.559 0.198 -0.849 -0.423     11.5 0.96
#> c[2]    -0.47915 -0.487 0.113 -0.594 -0.351     16.4 0.97
#> c[3]    -0.43951 -0.407 0.152 -0.576 -0.285      8.9 1.34
#> c[4]    -0.28090 -0.228 0.148 -0.485 -0.112     11.0 1.03
#> c[5]    -0.19411 -0.175 0.094 -0.304 -0.099     18.2 1.00
#> c[6]    -0.08573 -0.098 0.127 -0.257  0.074      3.0 1.44
#> c[7]     0.30059  0.311 0.126  0.173  0.465      2.9 1.52
#> c[8]     0.25602  0.282 0.164  0.022  0.462      9.5 1.27
#> c[9]     0.44350  0.507 0.158  0.237  0.607      2.8 1.55
#> c[10]    0.35560  0.247 0.248  0.065  0.667      2.1 2.12
#> A[1,1]   1.33305  1.460 0.272  0.935  1.614      3.2 1.37
#> A[2,1]   0.33478  0.379 0.204  0.071  0.547      2.8 1.50
#> A[3,1]   0.23170  0.271 0.261 -0.114  0.538      2.2 2.12
#> A[4,1]  -0.19306 -0.180 0.133 -0.357 -0.047      6.0 1.19
#> A[5,1]  -0.12407 -0.111 0.172 -0.352  0.112      6.9 1.14
#> A[6,1]   0.06257  0.060 0.132 -0.143  0.237     26.0 0.99
#> A[7,1]  -0.21693 -0.251 0.179 -0.411  0.024      7.0 1.27
#> A[8,1]  -0.26832 -0.317 0.157 -0.440 -0.046      9.0 1.04
#> A[9,1]   0.28778  0.297 0.179  0.025  0.515      2.3 1.84
#> A[10,1] -0.00023  0.037 0.162 -0.195  0.170     12.0 0.97
#> A[2,2]   0.70885  0.678 0.187  0.517  0.941      3.6 1.29
#> A[3,2]   0.38025  0.371 0.174  0.161  0.597      8.4 1.09
#> A[4,2]   0.15580  0.128 0.181 -0.050  0.423     14.7 1.01
#> A[5,2]   0.55081  0.626 0.378 -0.056  0.961      2.1 2.12
#> A[6,2]   0.49227  0.566 0.348  0.069  0.870      2.1 2.12
#> A[7,2]   0.68732  0.729 0.267  0.303  0.932      8.3 1.15
#> A[8,2]   0.46554  0.395 0.230  0.279  0.776      5.2 1.19
#> A[9,2]   0.59997  0.615 0.154  0.407  0.761      3.0 1.49
#> A[10,2]  0.70715  0.757 0.377  0.240  1.277      2.1 2.12
#> A[3,3]   0.65682  0.623 0.179  0.432  0.861      8.6 1.19
#> A[4,3]   0.54910  0.434 0.309  0.265  0.952     10.8 1.02
#> A[5,3]   0.17844  0.177 0.304 -0.205  0.529      2.1 2.12
#> A[6,3]  -0.03184  0.018 0.287 -0.543  0.356      3.0 1.46
#> A[7,3]   0.30015  0.276 0.135  0.125  0.488     13.0 1.01
#> A[8,3]   0.78520  0.838 0.286  0.442  1.069      2.1 1.97
#> A[9,3]   0.69462  0.740 0.249  0.427  1.000      5.4 1.19
#> A[10,3]  0.61241  0.659 0.167  0.425  0.783     11.5 0.96
#> 
#> ess_bulk is the bulk effective sample size; rhat is the potential
#> scale reduction factor on split chains (Rhat = 1 at convergence).
#> Use summary() for the full set of statistics (incl. ess_tail).
```
