# ---------------------------------------------------------------------------
# plot.spifa()
# ---------------------------------------------------------------------------

test_that("plot.spifa() returns a two-panel patchwork combining trace and density", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  gg <- plot(samples, select = "c")
  expect_s3_class(gg, "patchwork")
  expect_equal(length(gg$patches$plots), 1)
})

test_that("plot.spifa() shows the same parameters, row for row, in both panels", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  # select = "c" matches exactly 10 -- the nshow (10) default, no subsampling
  gg <- plot(samples, select = "c")
  trace_pars <- levels(droplevels(gg$patches$plots[[1]]$data$parameter))
  density_pars <- levels(droplevels(gg$data$parameter))
  expect_equal(trace_pars, density_pars)
  expect_length(trace_pars, 10)

  # select = "A" matches more than nshow -- both panels must still agree
  # (plot_trace()/plot_density() each apply their own independent random
  # subsample when called separately with the same select)
  gg_a <- plot(samples, select = "A", nshow = 6)
  trace_pars_a <- levels(droplevels(gg_a$patches$plots[[1]]$data$parameter))
  density_pars_a <- levels(droplevels(gg_a$data$parameter))
  expect_equal(trace_pars_a, density_pars_a)
  expect_length(trace_pars_a, 6)
})

test_that("plot.spifa() defaults to select = c('c', 'A') and drops restricted A entries", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(1, 6, 10), 1] <- 0
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors, ngp = 0,
    niter = 20, thin = 1,
    constraints = list(discrimination = A))

  gg <- plot(samples, nshow = NULL)
  pars <- as.character(gg$data$parameter)
  expect_false(any(grepl("^A\\[list\\(1,1\\)\\]$", pars)))
  expect_true(any(grepl("^c\\[", pars)))
  expect_true(any(grepl("^A\\[", pars)))
})

test_that("plot.spifa() nshow = NULL shows every matched parameter", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  gg <- plot(samples, select = "c", nshow = NULL)
  expect_length(levels(droplevels(gg$data$parameter)), 10)
})
