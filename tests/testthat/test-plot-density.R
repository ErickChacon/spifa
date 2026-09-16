# ---------------------------------------------------------------------------
# plot_density()
# ---------------------------------------------------------------------------

test_that("plot_density(): selects by group name or explicit indexed names", {
  data(ipixuna, package = "spifa")
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 20)

  gg_group <- plot_density(samples, select = "c", facet = TRUE)
  expect_s3_class(gg_group, "ggplot")
  expect_equal(length(unique(gg_group$data$parameter)), ncol(ipixuna$items))

  gg_explicit <- plot_density(samples, select = paste0("c[", 1:3, "]"), facet = FALSE)
  expect_equal(sort(as.character(unique(gg_explicit$data$parameter))),
    paste0("c[", 1:3, "]"))
})

test_that("plot_density(): orders parameters by index, not alphabetically", {
  data(ipixuna, package = "spifa")
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 20)

  # "c[10]" should sort after "c[9]", not between "c[1]" and "c[2]"
  gg <- plot_density(samples, select = "c", facet = TRUE, nshow = NULL)
  expect_equal(levels(gg$data$parameter), paste0("c[", 1:10, "]"))
})

test_that("plot_density(): drops restricted A parameters", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(1, 6, 10), 1] <- 0
  A[c(3, 5, 8), 2] <- 0
  A[c(2, 9), 3] <- 0

  samples <- spifa(
    items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 20,
    constraints = list(discrimination = A)
  )

  gg <- plot_density(samples, select = "A", facet = TRUE, nshow = NULL)
  expect_equal(length(unique(gg$data$parameter)), sum(A == 1))
})

test_that("plot_density(): nshow randomly subsamples parameters", {
  data(ipixuna, package = "spifa")
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 20)

  gg <- plot_density(samples, select = "c", facet = TRUE, nshow = 3)
  expect_equal(length(unique(gg$data$parameter)), 3)

  # nshow = NULL disables the cap, showing every matched parameter
  gg_all <- plot_density(samples, select = "c", facet = TRUE, nshow = NULL)
  expect_equal(length(unique(gg_all$data$parameter)), ncol(ipixuna$items))
})

test_that("plot_density(): facet respects ncol and facet_scales", {
  data(ipixuna, package = "spifa")
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 20)

  gg_default <- plot_density(samples, select = "c", facet = TRUE)
  expect_equal(gg_default$facet$params$ncol, 1)
  expect_true(gg_default$facet$params$free$x)

  gg_wide <- plot_density(samples, select = "c", facet = TRUE, ncol = 2)
  expect_equal(gg_wide$facet$params$ncol, 2)

  # facet_scales allowed
  gg_sharedx <- plot_density(samples, select = "c", facet = TRUE, facet_scales = "free_y")
  expect_false(gg_sharedx$facet$params$free$x)
  expect_true(gg_sharedx$facet$params$free$y)
})

test_that("plot_density(): plotmath-parses single and two-index parameter names", {
  data(ipixuna, package = "spifa")

  samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 1, niter = 20)

  # two-index "A[i,j]" becomes "A[list(i,j)]" for parsing and order
  # drop the second index
  gg_a <- plot_density(samples, select = "A[1,1]", facet = TRUE)
  expect_equal(as.character(unique(gg_a$data$parameter)), "A[list(1,1)]")

  # single-index "phi[1]" is left as-is (already valid, Greek-letter, plotmath)
  gg_phi <- plot_density(samples, select = "phi", facet = TRUE)
  expect_equal(as.character(unique(gg_phi$data$parameter)), "phi[1]")
})
