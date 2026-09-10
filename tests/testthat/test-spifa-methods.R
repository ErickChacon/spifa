# ---------------------------------------------------------------------------
# predict.spifa(): behaviour genuinely differs by model_type
# ---------------------------------------------------------------------------

test_that("predict.spifa(): eifa/cifa", {
 # returns the latent abilities of the observed subjects
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3

  # eifa
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5)
  pr <- predict(samples)
  expect_equal(pr, as_draws_array(as.list(samples)$Theta))
  expect_equal(predict(samples, newdata = sf::st_geometry(ipixuna)[1:2]), pr)

  # cifa
  A <- matrix(1, nitems, nfactors)
  A[c(1, 3, 5), 1] <- 0
  A[c(2, 6, 9), 2] <- 0
  A[c(4, 7, 10), 3] <- 0
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, constraints = list(discrimination = A))
  pr <- predict(samples)
  expect_equal(pr, as_draws_array(as.list(samples)$Theta))
  expect_equal(predict(samples, newdata = sf::st_geometry(ipixuna)[1:2]), pr)

  # burnin/thin subset the returned samples
  pr_burnin <- predict(samples, burnin = 2, thin = 2)
  expect_equal(pr_burnin,
    as_draws_array(as.list(samples)$Theta[c(3, 5), , drop = FALSE]))
})

test_that("predict.spifa(): spifa", {
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
  expect_equal(pr, as_draws_array(as.list(samples)$Theta))

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
})

test_that("predict.spifa(): cifa_pred", {
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
  expect_equal(pr_none, as_draws_array(as.list(samples)$Theta))

  # incorrect newdata: missing predictors
  expect_error(predict(samples, newdata = data.frame(other_col = c(1, 2, 3))),
    "missing the predictor column")

  # adequate newdata
  newdata <- data.frame(wealth = c(0.1, -0.2, 0.3))
  pr <- predict(samples, newdata = newdata)
  expect_s3_class(pr, "draws_array")
  expect_equal(ndraws(pr), 5)
  expect_equal(nvariables(pr), 3 * nfactors)
})

test_that("predict.spifa(): spifa_pred", {
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
  expect_equal(pr_none, as_draws_array(as.list(samples)$Theta))

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
})

# ---------------------------------------------------------------------------
# dic/summary/as_tibble/as.list/print: model_type-agnostic, tested once
# ---------------------------------------------------------------------------

test_that("dic.spifa() computes deviance information criterion components", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(3, 9), 1] <- 0
  A[c(1, 5, 8), 2] <- 0
  A[c(2, 6, 7, 10), 3] <- 0
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1,
    constraints = list(discrimination = A))

  d <- dic(samples)
  expect_true(all(c("average_of_deviance", "n_effec_params", "dic") %in% names(d)))
})

test_that("summary.spifa() returns posterior summaries with burnin/thin/select", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(1, 3, 5), 1] <- 0
  A[c(2, 6, 9), 2] <- 0
  A[c(4, 7, 10), 3] <- 0
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 10, thin = 1,
    constraints = list(discrimination = A))

  smry <- summary(samples, select = "c")
  expect_true(all(c("variable", "q2.5", "median", "q97.5") %in% names(smry)))
  expect_equal(nrow(smry), ncol(ipixuna_flat$items))

  smry_burnin <- summary(samples, burnin = 5, select = "c")
  c_block <- as.list(samples)$c
  expect_equal(smry_burnin$mean, unname(colMeans(c_block[6:10, , drop = FALSE])))
})

test_that("as_tibble.spifa() converts to a wide tibble with burnin/thin/select", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(2, 7), 1] <- 0
  A[c(1, 4, 9, 10), 2] <- 0
  A[c(3, 6), 3] <- 0
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 10, thin = 1,
    constraints = list(discrimination = A))

  samples_tib <- as_tibble(samples)
  expect_false(inherits(samples_tib, "spifa"))
  expect_equal(nrow(samples_tib), 10)

  samples_tib_sub <- as_tibble(samples, burnin = 5, thin = 2, select = "c")
  expect_equal(nrow(samples_tib_sub), length(seq(6, 10, 2)))
  expect_true(all(grepl("^c\\[", names(samples_tib_sub))))
})

test_that("as.list.spifa() splits samples back into block-shaped matrices", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(1, 6, 10), 1] <- 0
  A[c(3, 5, 8), 2] <- 0
  A[c(2, 9), 3] <- 0
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1,
    constraints = list(discrimination = A))

  samples_list <- as.list(samples)
  expect_equal(nrow(samples_list$c), 5)
  expect_equal(ncol(samples_list$c), ncol(ipixuna_flat$items))
})

test_that("print.spifa() prints without error", {
  data(ipixuna, package = "spifa")
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0,
    niter = 5, thin = 1)

  expect_output(print(samples))
})
