# Tests for the S3 methods on a fitted spifa object (R/spifa-methods.R).
# predict.spifa() genuinely branches on model_type, so it gets one test per
# branch; dic/summary/as_tibble/as.list/print don't inspect model_type at
# all, so each is tested once against whichever fitted object is cheapest.
# See test-spifa.R for spifa()'s own fitting behaviour.

# ---------------------------------------------------------------------------
# predict.spifa(): behaviour genuinely differs by model_type
# ---------------------------------------------------------------------------

test_that("predict.spifa() returns the means-only message for eifa/cifa", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  expect_equal(predict(samples), "here is only the means of the latent abilities")
  # newcoords/newdata are meaningless for eifa/cifa -- still short-circuits
  expect_equal(predict(samples, newcoords = matrix(1:4, 2, 2)),
    "here is only the means of the latent abilities")
})

test_that("predict.spifa() predicts the spatial process at newcoords for spifa", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)

  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, mgp = diag(nfactors), resid_sd = parameters$resid_params$sd),
    priors = list(
      mgp_sd = list(initial = 0.6, mean = 0.6, sd = 0.4),
      mgp_range = list(initial = 200, mean = 200, sd = 0.4)))

  newcoords <- sf::st_coordinates(ipixuna$geometry)[1:5, , drop = FALSE]
  pr <- predict(samples, newcoords = newcoords)
  expect_type(pr, "list")
  expect_equal(dim(pr$theta), c(5, 5 * nfactors))
})

test_that("predict.spifa() predicts new subjects from newdata for cifa_pred", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ wealth, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  # one column: formula's right-hand side (~ wealth) has no intercept, just
  # the single "wealth" predictor
  newdata <- matrix(c(0.1, -0.2, 0.3), 3, 1)
  pr <- predict(samples, newdata = newdata)
  expect_type(pr, "list")
  expect_equal(dim(pr$theta), c(5, 3 * nfactors))
})

test_that("predict.spifa() falls back to reference-level predictors for spifa_pred without newdata", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)

  samples <- spifa(items ~ wealth, data = ipixuna, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, mgp = diag(nfactors), resid_sd = parameters$resid_params$sd),
    priors = list(
      mgp_sd = list(initial = 0.6, mean = 0.6, sd = 0.4),
      mgp_range = list(initial = 200, mean = 200, sd = 0.4)))

  newcoords <- sf::st_coordinates(ipixuna$geometry)[1:5, , drop = FALSE]
  pr <- predict(samples, newcoords = newcoords)
  expect_type(pr, "list")
  expect_equal(dim(pr$theta), c(5, 5 * nfactors))
})

# ---------------------------------------------------------------------------
# dic/summary/as_tibble/as.list/print: model_type-agnostic, tested once
# ---------------------------------------------------------------------------

test_that("dic.spifa() computes deviance information criterion components", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  d <- dic(samples)
  expect_true(all(c("average_of_deviance", "n_effec_params", "dic") %in% names(d)))
})

test_that("summary.spifa() returns posterior summaries with burnin/thin/select", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 10, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  smry <- summary(samples, select = "c")
  expect_true(all(c("variable", "q2.5", "median", "q97.5") %in% names(smry)))
  expect_equal(nrow(smry), ncol(ipixuna_flat$items))

  smry_burnin <- summary(samples, burnin = 5, select = "c")
  c_block <- as.list(samples)$c
  expect_equal(smry_burnin$mean, unname(colMeans(c_block[6:10, , drop = FALSE])))
})

test_that("as_tibble.spifa() converts to a wide tibble with burnin/thin/select", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 10, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  samples_tib <- as_tibble(samples)
  expect_false(inherits(samples_tib, "spifa"))
  expect_equal(nrow(samples_tib), 10)

  samples_tib_sub <- as_tibble(samples, burnin = 5, thin = 2, select = "c")
  expect_equal(nrow(samples_tib_sub), length(seq(6, 10, 2)))
  expect_true(all(grepl("^c\\[", names(samples_tib_sub))))
})

test_that("as.list.spifa() splits samples back into block-shaped matrices", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  samples_list <- as.list(samples)
  expect_equal(nrow(samples_list$c), 5)
  expect_equal(ncol(samples_list$c), ncol(ipixuna_flat$items))
})

test_that("print.spifa() prints without error", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  nfactors <- ncol(parameters$discrimination)
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1, standardize = FALSE)

  expect_output(print(samples))
})
