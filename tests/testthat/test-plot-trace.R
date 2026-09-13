# ---------------------------------------------------------------------------
# plot_trace()
# ---------------------------------------------------------------------------

test_that("plot_trace() selects parameters by block name or explicit indexed names", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  gg_block <- plot_trace(samples, select = "c", facet = TRUE)
  expect_s3_class(gg_block, "ggplot")
  expect_equal(length(unique(gg_block$data$parameter)), ncol(ipixuna_flat$items))

  gg_explicit <- plot_trace(samples, select = paste0("c[", 1:3, "]"), facet = FALSE)
  expect_equal(sort(as.character(unique(gg_explicit$data$parameter))), paste0("c[", 1:3, "]"))
})

test_that("plot_trace() orders parameters naturally, not alphabetically", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  # "c[10]" should sort after "c[9]", not between "c[1]" and "c[2]"
  gg <- plot_trace(samples, select = "c", facet = TRUE, nshow = NULL)
  expect_equal(levels(gg$data$parameter), paste0("c[", 1:10, "]"))
})

test_that("plot_trace() drops structurally-restricted A parameters via drop_restricted()", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(1, 6, 10), 1] <- 0
  A[c(3, 5, 8), 2] <- 0
  A[c(2, 9), 3] <- 0
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors, ngp = 0,
    niter = 20, thin = 1, constraints = list(discrimination = A))

  gg <- plot_trace(samples, select = "A", facet = TRUE, nshow = NULL)
  expect_equal(length(unique(gg$data$parameter)), sum(A == 1))
})

test_that("plot_trace() nshow randomly subsamples when select matches more parameters", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  gg <- plot_trace(samples, select = "c", facet = TRUE, nshow = 3)
  expect_equal(length(unique(gg$data$parameter)), 3)

  # nshow = NULL disables the cap, showing every matched parameter
  gg_all <- plot_trace(samples, select = "c", facet = TRUE, nshow = NULL)
  expect_equal(length(unique(gg_all$data$parameter)), ncol(ipixuna_flat$items))
})

test_that("plot_trace() facet format respects ncol", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  gg_default <- plot_trace(samples, select = "c", facet = TRUE)
  expect_equal(gg_default$facet$params$ncol, 1)

  gg_wide <- plot_trace(samples, select = "c", facet = TRUE, ncol = 2)
  expect_equal(gg_wide$facet$params$ncol, 2)
})

test_that("plot_trace() overlay legend auto-hides past 10 series and can be overridden", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors, ngp = 0,
    niter = 20, thin = 1, constraints = list(discrimination = A))

  # c has 10 series (<= 10): legend shown by default
  gg_c <- plot_trace(samples, select = "c", facet = FALSE)
  expect_equal(gg_c$theme$legend.position, "bottom")

  # A (unrestricted here) has 30 series (> 10): legend hidden by default
  gg_a <- plot_trace(samples, select = "A", facet = FALSE, nshow = NULL)
  expect_equal(gg_a$theme$legend.position, "none")

  # explicit override in either direction
  gg_a_forced <- plot_trace(samples, select = "A", facet = FALSE,
    nshow = NULL, legend = "bottom")
  expect_equal(gg_a_forced$theme$legend.position, "bottom")

  gg_c_hidden <- plot_trace(samples, select = "c", facet = FALSE, legend = FALSE)
  expect_equal(gg_c_hidden$theme$legend.position, "none")
})

test_that("plot_trace() plotmath-parses two-index and single-index parameter names", {
  data(ipixuna, package = "spifa")

  # ngp = 1 needs spatial coordinates, so this one keeps the sf geometry
  # (unlike the other tests here, which fit a plain, non-spatial model)
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 1,
    niter = 20, thin = 1)

  # two-index "A[i,j]" becomes "A[list(i,j)]" so plotmath doesn't silently
  # drop the second index
  gg_a <- plot_trace(samples, select = "A[1,1]", facet = TRUE)
  expect_equal(as.character(unique(gg_a$data$parameter)), "A[list(1,1)]")

  # single-index "phi[1]" is left as-is (already valid, Greek-letter, plotmath)
  gg_phi <- plot_trace(samples, select = "phi", facet = TRUE)
  expect_equal(as.character(unique(gg_phi$data$parameter)), "phi[1]")
})
