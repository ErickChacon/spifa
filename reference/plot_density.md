# Density Plot of spifa Posterior Samples

Draws posterior density curves directly from a fitted `spifa` model, in
one of two formats: overlaid (`facet = FALSE`, the default; stacked,
slightly-overlapping ridgeline plot via
[`geom_density_ridges`](https://wilkelab.org/ggridges/reference/geom_density_ridges.html)
– the more insightful default for comparing many parameters' shapes and
locations at a glance) or faceted (`facet = TRUE`; one panel per
parameter). Shows the density shape only – for credible intervals and
point estimates, see
[`plot_interval`](https://ErickChacon.github.io/spifa/reference/plot_interval.md).
The density itself is computed directly by
[`geom_density`](https://ggplot2.tidyverse.org/reference/geom_density.html)/[`geom_density_ridges`](https://wilkelab.org/ggridges/reference/geom_density_ridges.html)
from the raw draws, reshaped via
[`as_draws_df`](https://mc-stan.org/posterior/reference/draws_df.html).

## Usage

``` r
plot_density(
  x,
  select,
  facet = FALSE,
  burnin = 0,
  thin = 1,
  nshow = 10,
  ncol = 1,
  facet_scales = "free",
  scale = 1.2,
  ...
)
```

## Arguments

- x:

  A fitted `spifa` model.

- select:

  Parameters to plot, as in
  [`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md).

- facet:

  Logical; if `TRUE`, draw one panel per parameter; if `FALSE`
  (default), overlay every parameter as a ridgeline plot. See
  Description.

- burnin:

  Number of initial iterations to discard.

- thin:

  Thinning interval applied after `burnin`.

- nshow:

  As in
  [`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md).

- ncol:

  Number of columns in the facet grid (`facet = TRUE` only); defaults to
  a single column.

- facet_scales:

  The `scales` argument of
  [`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  (`facet = TRUE` only): one of `"free"` (default; each panel gets its
  own x/y-axis), `"free_y"` (one shared x-axis, line shown only on the
  bottom panel of each column, matching
  [`plot_trace`](https://ErickChacon.github.io/spifa/reference/plot_trace.md)'s
  facet – useful for comparing parameters on a similar scale),
  `"free_x"`, or `"fixed"`. Parameters with very different value ranges
  can render compressed/hard to read with a shared x-axis.

- scale:

  Amount of vertical overlap between ridges (`facet = FALSE` only),
  passed to
  [`geom_density_ridges`](https://wilkelab.org/ggridges/reference/geom_density_ridges.html).
  Defaults to `1.2`.

- ...:

  Currently unused.

## Value

A `ggplot` object.

## Author

Erick A. Chacón-Montalván

## Examples

``` r
# \donttest{
data(ipixuna)
samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 1000)

plot_density(samples, select = "c", facet = FALSE)
#> Picking joint bandwidth of 0.0509

plot_density(samples, select = "c", facet = TRUE)


# more than nshow (10) parameters: a random subsample is shown
plot_density(samples, select = "A", facet = FALSE, nshow = 6)
#> Picking joint bandwidth of 0.0973

# }
```
