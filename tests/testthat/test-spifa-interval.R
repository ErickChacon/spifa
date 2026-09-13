# ---------------------------------------------------------------------------
# plot_interval()
# ---------------------------------------------------------------------------

test_that("plot_interval() selects parameters by block name or explicit indexed names", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  gg_block <- plot_interval(samples, select = "c")
  expect_s3_class(gg_block, "ggplot")
  expect_equal(length(unique(gg_block$data$parameter)), ncol(ipixuna_flat$items))

  gg_explicit <- plot_interval(samples, select = paste0("c[", 1:3, "]"))
  expect_equal(sort(as.character(unique(gg_explicit$data$parameter))), paste0("c[", 1:3, "]"))
})

test_that("plot_interval() orders parameters naturally, not alphabetically", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  # "c[10]" should sort after "c[9]", not between "c[1]" and "c[2]"
  gg <- plot_interval(samples, select = "c")
  expect_equal(levels(gg$data$parameter), paste0("c[", 1:10, "]"))
})

test_that("plot_interval() drops structurally-restricted A parameters via drop_restricted()", {
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

  gg <- plot_interval(samples, select = "A")
  expect_equal(length(unique(gg$data$parameter)), sum(A == 1))
})

test_that("plot_interval() nshow randomly subsamples when select matches more parameters", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  gg <- plot_interval(samples, select = "c", nshow = 3)
  expect_equal(length(unique(gg$data$parameter)), 3)

  # NULL (the default) shows every matched parameter, unlike
  # plot_trace()/plot_density(), which cap at 10 by default
  gg_all <- plot_interval(samples, select = "c")
  expect_equal(length(unique(gg_all$data$parameter)), ncol(ipixuna_flat$items))
})

test_that("plot_interval() horizontal swaps which axis carries the parameters", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  # default (FALSE): parameters on x, values on y -- the SPIFA paper's own
  # ci_intervals() convention
  gg_vertical <- plot_interval(samples, select = "c")
  expect_null(gg_vertical$labels$x)
  expect_equal(gg_vertical$labels$y, "Value")

  # TRUE: axes swapped (bayesplot::mcmc_intervals()'s own layout)
  gg_horizontal <- plot_interval(samples, select = "c", horizontal = TRUE)
  expect_equal(gg_horizontal$labels$x, "Value")
  expect_null(gg_horizontal$labels$y)
})

test_that("plot_interval() sort reorders parameters by point estimate", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  gg <- plot_interval(samples, select = "c", sort = TRUE)
  point_est <- gg$data$m[match(levels(gg$data$parameter), gg$data$parameter)]
  expect_equal(point_est, sort(point_est))
})

test_that("plot_interval() plotmath-parses two-index and single-index parameter names", {
  data(ipixuna, package = "spifa")

  # ngp = 1 needs spatial coordinates, so this one keeps the sf geometry
  # (unlike the other tests here, which fit a plain, non-spatial model)
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 1,
    niter = 20, thin = 1)

  # two-index "A[i,j]" becomes "A[list(i,j)]" so plotmath doesn't silently
  # drop the second index
  gg_a <- plot_interval(samples, select = "A[1,1]")
  expect_equal(as.character(unique(gg_a$data$parameter)), "A[list(1,1)]")

  # single-index "phi[1]" is left as-is (already valid, Greek-letter, plotmath)
  gg_phi <- plot_interval(samples, select = "phi")
  expect_equal(as.character(unique(gg_phi$data$parameter)), "phi[1]")
})

test_that("plot_interval() reference accepts a matrix matching a block's own shape", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  nitems <- ncol(ipixuna_flat$items)
  nfactors <- 3
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors, ngp = 0,
    niter = 20, thin = 1)

  ref_mat <- matrix(seq_len(nitems * nfactors), nitems, nfactors)
  gg <- plot_interval(samples, select = "A", reference = ref_mat)

  expect_true("reference" %in% names(gg$data))
  # data.frame is keyed by raw "A[i,j]" names before plotmath parsing, so
  # check the join picked up the right matrix entry per parameter
  raw <- gsub("list\\((.+),(.+)\\)", "\\1,\\2", as.character(gg$data$parameter))
  raw <- gsub("^A\\[|\\]$", "", raw)
  idx <- do.call(rbind, strsplit(raw, ","))
  expected <- ref_mat[cbind(as.integer(idx[, 1]), as.integer(idx[, 2]))]
  expect_equal(gg$data$reference, expected)
})

test_that("plot_interval() reference accepts a plain vector matching a block's own shape", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  nitems <- ncol(ipixuna_flat$items)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  ref <- seq_len(nitems) / 2
  gg <- plot_interval(samples, select = "c", reference = ref)

  by_param <- setNames(gg$data$reference, as.character(gg$data$parameter))
  expect_equal(unname(by_param[paste0("c[", seq_len(nitems), "]")]), ref)
})

test_that("plot_interval() reference requires select to be a single block name", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  expect_error(
    plot_interval(samples, select = "c[1]", reference = 0.5),
    "single block")
})

test_that("plot_interval() reference interacts correctly with sort and horizontal", {
  data(ipixuna, package = "spifa")
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  nitems <- ncol(ipixuna_flat$items)
  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = 3, ngp = 0,
    niter = 20, thin = 1)

  ref <- seq_len(nitems)
  gg <- plot_interval(samples, select = "c", reference = ref, sort = TRUE,
    horizontal = TRUE)

  expect_s3_class(gg, "ggplot")
  expect_equal(nrow(gg$data), nitems)
  # reference values still line up with their own parameter after reordering
  by_param <- setNames(gg$data$reference, as.character(gg$data$parameter))
  expect_equal(unname(by_param[paste0("c[", seq_len(nitems), "]")]), ref)
})
