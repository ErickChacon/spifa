test_that("spifa() fits a confirmatory (non-spatial) model end-to-end", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)

  samples <- spifa(
    items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  expect_s3_class(samples, "spifa")
  expect_s3_class(samples, "draws_array")
  samples_list <- as_list(samples)
  # cifa: no predictors, no spatial structure -- mgp_sd/mgp_phi/betas were
  # never sampled and are dropped from the output (see R/spifa.R).
  expect_setequal(names(samples_list),
    c("c", "a", "theta", "z", "corr_chol", "corr"))
  expect_equal(nrow(samples_list$c), 5)
  expect_equal(ncol(samples_list$c), 10)

  samples_tib <- as_tibble(samples)
  expect_false(inherits(samples_tib, "spifa"))
  expect_equal(nrow(samples_tib), 5)

  smry <- summary(samples, select = "c")
  expect_true(all(c("variable", "q2.5", "median", "q97.5") %in% names(smry)))

  d <- dic(samples)
  expect_true(all(c("average_of_deviance", "n_effec_params", "dic") %in% names(d)))
})

test_that("spifa() fits a spatial model with predictors and predicts at new locations", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)

  samples <- spifa(
    items ~ wealth, data = ipixuna, nfactors = nfactors, niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, mgp = diag(nfactors), resid_sd = parameters$resid_params$sd),
    priors = list(
      mgp_sd = list(initial = 0.6, mean = 0.6, sd = 0.4),
      mgp_range = list(initial = 200, mean = 200, sd = 0.4)))

  expect_s3_class(samples, "spifa")
  # spifa_pred: has predictors and spatial structure -- every block was
  # actually sampled and none should be dropped.
  expect_setequal(unique(.spifa_var_blocks(dimnames(samples)[[3]])),
    c("c", "a", "theta", "z", "corr_chol", "corr", "mgp_sd", "mgp_phi", "betas"))

  newcoords <- sf::st_coordinates(ipixuna$geometry)[1:5, , drop = FALSE]
  pr <- predict(samples, newcoords = newcoords)
  expect_type(pr, "list")
  expect_equal(dim(pr$theta), c(5, 5 * nfactors))
})
