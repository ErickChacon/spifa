# ---------------------------------------------------------------------------
# dic.spifa()/predict.spifa(): both related to C++ functions
# ---------------------------------------------------------------------------

# dic() check
expect_dic_ok <- function (samples) {
  d <- dic(samples)
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 1)
  expect_true(all(c("mean_deviance", "p_eff", "dic") %in% names(d)))
  invisible(d)
}

test_that("predict/dic: eifa/cifa", {
 # returns the latent abilities of the observed subjects
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3

  # eifa
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5)
  pr <- predict(samples)
  expect_equal(pr, as_draws_array(subset_draws(samples, variable = "Theta")), ignore_attr = "class")
  expect_equal(predict(samples, newdata = sf::st_geometry(ipixuna)[1:2]), pr)
  expect_dic_ok(samples)

  # cifa
  A <- matrix(1, nitems, nfactors)
  A[c(1, 3, 5), 1] <- 0
  A[c(2, 6, 9), 2] <- 0
  A[c(4, 7, 10), 3] <- 0
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, constraints = list(discrimination = A))
  pr <- predict(samples)
  expect_equal(pr, as_draws_array(subset_draws(samples, variable = "Theta")), ignore_attr = "class")
  expect_equal(predict(samples, newdata = sf::st_geometry(ipixuna)[1:2]), pr)
  expect_dic_ok(samples)

  # burnin/thin subset the returned samples
  pr_burnin <- predict(samples, burnin = 2, thin = 2)
  expect_equal(pr_burnin,
    as_draws_array(subset_draws(samples, variable = "Theta", iteration = c(3, 5))),
    ignore_attr = "class")
})

test_that("predict/dic: spifa", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0

  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
    niter = 5, constraints = list(discrimination = A, loading = diag(nfactors)),
    priors = list(
      loading = list(initial = 0.6, mean = 0.6, sd = 0.4),
      range = list(initial = 200, mean = 200, sd = 0.4))
  )

  # no newdata: latent abilities at observed locations
  pr <- predict(samples)
  expect_equal(pr, as_draws_array(subset_draws(samples, variable = "Theta")), ignore_attr = "class")

  # newdata (sfc): latent abilities at new locations
  newcoords <- sf::st_make_grid(ipixuna, n = c(5, 1), what = "centers")
  # marginal spatial prediction
  pr <- predict(samples, newdata = newcoords)
  expect_s3_class(pr, "draws_array")
  expect_equal(ndraws(pr), 5)
  expect_equal(nvariables(pr), 5 * nfactors)
  # joint spatial prediction
  pr_joint <- predict(samples, newdata = newcoords, joint = TRUE)
  expect_equal(dim(pr_joint), dim(pr))

  expect_dic_ok(samples)
})

test_that("predict/dic: cifa_pred", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(2, 7), 1] <- 0
  A[c(1, 4, 9, 10), 2] <- 0
  A[c(3, 6), 3] <- 0

  samples <- spifa(items ~ wealth, data = ipixuna, nfactors = nfactors,
    ngp = 0, niter = 5, thin = 1, constraints = list(discrimination = A))

  # no newdata: latent abilities at observed subjects
  pr_none <- predict(samples)
  expect_equal(pr_none, as_draws_array(subset_draws(samples, variable = "Theta")), ignore_attr = "class")

  # incorrect newdata: missing predictors
  expect_error(predict(samples, newdata = data.frame(other_col = c(1, 2, 3))),
    "missing the predictor column")

  # adequate newdata
  newdata <- data.frame(wealth = c(0.1, -0.2, 0.3))
  pr <- predict(samples, newdata = newdata)
  expect_s3_class(pr, "draws_array")
  expect_equal(ndraws(pr), 5)
  expect_equal(nvariables(pr), 3 * nfactors)

  expect_dic_ok(samples)
})

test_that("predict/dic: spifa_pred", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0
  ngp <- 2
  loading <- matrix(c(1, 1, 0, 0, 0, 1), nfactors, ngp)

  samples <- spifa(items ~ wealth, data = ipixuna, nfactors = nfactors, ngp = ngp,
    niter = 5, constraints = list(discrimination = A, loading = loading),
    priors = list(
      loading = list(initial = 0.6, mean = 0.6, sd = 0.4),
      range = list(initial = 200, mean = 200, sd = 0.4))
  )

  # no newdata: latent abilities at observed locations
  pr_none <- predict(samples)
  expect_equal(pr_none, as_draws_array(subset_draws(samples, variable = "Theta")), ignore_attr = "class")

  # incorrect newdata: missing predictors
  newcoords <- sf::st_make_grid(ipixuna, n = c(5, 1), what = "centers")
  expect_error(predict(samples, newdata = newcoords), "wealth")
  expect_error(predict(samples, newdata = sf::st_sf(geometry = newcoords)), "wealth")
  expect_error(predict(samples, newdata = sf::st_sf(other_col = 1:5, geometry = newcoords)), "wealth")

  # incorrect newdata: a plain data.frame (no geometry) for a spatial model
  expect_error(predict(samples, newdata = data.frame(wealth = c(0.1, -0.2, 0.3))),
    "sf/sfc")

  # adequate newdata: supplying the predictor explicitly
  newdata <- sf::st_sf(wealth = rep(0, 5), geometry = newcoords)
  pr <- predict(samples, newdata = newdata)
  expect_s3_class(pr, "draws_array")
  expect_equal(ndraws(pr), 5)
  expect_equal(nvariables(pr), 5 * nfactors)

  expect_dic_ok(samples)
})

test_that("predict()/dic(): burnin/thin", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(3, 9), 1] <- 0
  A[c(1, 5, 8), 2] <- 0
  A[c(2, 6, 7, 10), 3] <- 0

  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 10, thin = 1, constraints = list(discrimination = A))
  d <- expect_dic_ok(samples)

  # dic filtering samples
  d_burnin <- dic(samples, burnin = 5, thin = 2)
  expect_false(isTRUE(all.equal(d, d_burnin)))

  # expected errors: burnin >= niter
  expect_error(predict(samples, burnin = 10), "iterations")
  expect_error(dic(samples, burnin = 15), "iterations")
})
