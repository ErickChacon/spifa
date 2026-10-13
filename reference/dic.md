# Deviance Information Criterion

Generic function for the Deviance Information Criterion (DIC), a
Bayesian measure of model fit that penalises complexity. See
[`dic.spifa`](https://ErickChacon.github.io/spifa/reference/dic.spifa.md)
for the `spifa` method.

## Usage

``` r
dic(x, ...)
```

## Arguments

- x:

  A fitted model object.

- ...:

  Further arguments passed to methods.

## Value

The DIC for `x`; see the specific method (e.g.
[`dic.spifa`](https://ErickChacon.github.io/spifa/reference/dic.spifa.md))
for its exact return value.
