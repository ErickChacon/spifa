test_that("spifa() fits a confirmatory (non-spatial) model end-to-end", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna_wide, "parameters")
  L_a <- (parameters$discrimination != 0) * 1

  samples <- spifa(
    responses = `Item 1`:`Item 10`, data = ipixuna_wide, nfactors = 2,
    niter = 5, thin = 1, standardize = FALSE,
    constrains = list(A = L_a, W = diag(2), V_sd = rep(0.5, 2)))

  expect_s3_class(samples, "spifa.list")
  expect_setequal(names(samples),
    c("c", "a", "theta", "z", "corr_chol", "corr", "mgp_sd", "mgp_phi", "betas"))
  expect_equal(nrow(samples$c), 5)
  expect_equal(ncol(samples$c), 10)

  samples_tib <- as_tibble(samples)
  expect_s3_class(samples_tib, "spifa")
  expect_equal(nrow(samples_tib), 5)

  smry <- summary(as_tibble(samples, select = "c"))
  expect_true(all(c("Parameters", "2.5%", "50%", "97.5%") %in% names(smry)))

  d <- dic(samples_tib)
  expect_true(all(c("average_of_deviance", "n_effec_params", "dic") %in% names(d)))
})

test_that("spifa() fits a spatial model with predictors and predicts at new locations", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna_wide, "parameters")
  L_a <- (parameters$discrimination != 0) * 1

  samples <- spifa(
    responses = `Item 1`:`Item 10`, pred_formula = ~ x1, coords = coords,
    data = ipixuna_wide, nfactors = 2, niter = 5, thin = 1, standardize = FALSE,
    constrains = list(A = L_a, W = diag(2), V_sd = rep(0.4^0.5, 2)),
    sigmas_gp_opt = list(initial = 0.6, prior_mean = 0.6, prior_sd = 0.4),
    phi_gp_opt = list(initial = 200, prior_mean = 200, prior_sd = 0.4))

  expect_s3_class(samples, "spifa.list")

  samples_tib <- as_tibble(samples)
  newcoords <- sf::st_coordinates(ipixuna_wide$coords)[1:5, , drop = FALSE]
  pr <- predict(samples_tib, newcoords = newcoords)
  expect_type(pr, "list")
  expect_equal(dim(pr$theta), c(5, 5 * 2))
})
