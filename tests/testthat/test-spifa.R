# Tests for spifa() itself: the model-type detection matrix (type-dependent
# behaviour) followed by sampler mechanics and argument validation (shared,
# model-type-agnostic behaviour). Method tests (predict/dic/summary/
# as_tibble/as.list/print) live in test-spifa-methods.R instead -- see
# dev/design/scope.md for the model_type definitions.

# ---------------------------------------------------------------------------
# Per-model-type: block presence, constraints/priors application, standardize
# ---------------------------------------------------------------------------

test_that("spifa() fits eifa (exploratory, no constraints) correctly", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(
    items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE)

  expect_equal(attr(samples, "spifa_args")$model_type, "eifa")
  expect_setequal(names(as.list(samples)),
    c("c", "a", "theta", "z", "corr_chol", "corr"))

  # standardize is a documented no-op for eifa/cifa (no GP/predictor scale
  # non-identifiability to normalize against)
  set.seed(42)
  samples_std <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, standardize = TRUE)
  set.seed(42)
  samples_raw <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE)
  strip_spifa_args <- function (x) { attr(x, "spifa_args") <- NULL; unclass(x) }
  expect_equal(strip_spifa_args(as.list(samples_std)), strip_spifa_args(as.list(samples_raw)))
})

test_that("spifa() fits cifa (confirmatory, no predictors/spatial) correctly", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(
    items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  expect_equal(attr(samples, "spifa_args")$model_type, "cifa")

  # constraints$discrimination mask respected: entries fixed to 0 stay
  # exactly 0 in every stored draw (parsed off the "A[i,j]" column names
  # rather than assumed, since the flatten order isn't a public contract)
  a_samples <- as.list(samples)$a
  idx <- as.integer(sub("A\\[(\\d+),(\\d+)\\]", "\\1", colnames(a_samples)))
  jdx <- as.integer(sub("A\\[(\\d+),(\\d+)\\]", "\\2", colnames(a_samples)))
  zero_mask <- L_a[cbind(idx, jdx)] == 0
  expect_true(all(a_samples[, zero_mask] == 0))
})

test_that("spifa() fits cifa_pred (confirmatory with predictors) correctly", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  effect_initial <- matrix(parameters$effect, nrow = 1)

  samples <- spifa(
    items ~ wealth, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd),
    priors = list(effect = list(initial = effect_initial, mean = effect_initial)))

  expect_equal(attr(samples, "spifa_args")$model_type, "cifa_pred")
  blocks <- as.list(samples)
  expect_setequal(names(blocks), c("c", "a", "theta", "z", "corr_chol", "corr", "betas"))
  expect_equal(dim(blocks$betas), c(5, ncol(effect_initial)))

  # standardize actually rescales theta to ~unit variance per factor (coarse
  # check: pool all iterations/observations for a factor and look at the sd
  # of the pooled vector, not an exact per-column target)
  set.seed(1)
  samples_std <- spifa(
    items ~ wealth, data = ipixuna_flat, nfactors = nfactors,
    niter = 200, thin = 1, standardize = TRUE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))
  theta <- as.list(samples_std)$theta
  factor_idx <- as.integer(sub(".*,(\\d+)\\]$", "\\1", colnames(theta)))
  for (k in seq_len(nfactors)) {
    expect_true(abs(sd(as.numeric(theta[, factor_idx == k])) - 1) < 0.5)
  }
})

test_that("spifa() fits spifa (spatial, no predictors) correctly", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)

  samples <- spifa(
    items ~ 1, data = ipixuna, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, mgp = diag(nfactors), resid_sd = parameters$resid_params$sd),
    priors = list(
      mgp_sd = list(initial = parameters$mgp_params$sd, mean = parameters$mgp_params$sd, sd = 0.4),
      mgp_range = list(initial = parameters$mgp_params$phi, mean = parameters$mgp_params$phi, sd = 0.3)))

  expect_equal(attr(samples, "spifa_args")$model_type, "spifa")
  blocks <- as.list(samples)
  expect_setequal(names(blocks),
    c("c", "a", "theta", "z", "corr_chol", "corr", "mgp_sd", "mgp_phi"))
  # constraints$mgp = diag(nfactors): one independent GP per factor, so
  # exactly nfactors mgp_sd entries are estimated (not the full nfactors x
  # ngp product)
  expect_equal(ncol(blocks$mgp_sd), nfactors)
})

test_that("spifa() fits spifa_pred (spatial with predictors) correctly", {
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

  expect_equal(attr(samples, "spifa_args")$model_type, "spifa_pred")
  expect_setequal(names(as.list(samples)),
    c("c", "a", "theta", "z", "corr_chol", "corr", "mgp_sd", "mgp_phi", "betas"))
})

test_that("spifa() with ngp = 0 opts an sf dataset out of the spatial model", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)

  # ipixuna is an sf object, but ngp = 0 must force a non-spatial fit anyway
  samples <- spifa(
    items ~ wealth, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1, standardize = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  expect_equal(attr(samples, "spifa_args")$model_type, "cifa_pred")
  expect_false(any(c("mgp_sd", "mgp_phi") %in% names(as.list(samples))))
})

# ---------------------------------------------------------------------------
# Shared, model-type-agnostic: sampler mechanics and argument validation
# ---------------------------------------------------------------------------

test_that("spifa() burnin discards iterations without storing them", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)
  constraints <- list(discrimination = L_a, resid_sd = parameters$resid_params$sd)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, burnin = 0, standardize = FALSE, constraints = constraints)
  expect_equal(dim(samples)[1], 5)

  samples_burnin <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, burnin = 10, standardize = FALSE, constraints = constraints)
  expect_equal(dim(samples_burnin)[1], 5)

  # non-divisible edge case: ceiling(niter / thin) draws are stored (this
  # combination previously under-allocated, see src/ifa.cpp's nsave fix)
  samples_edge <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 7, thin = 3, burnin = 2, standardize = FALSE, constraints = constraints)
  expect_equal(dim(samples_edge)[1], ceiling(7 / 3))
})

test_that("spifa() execute = FALSE returns spifa_args without sampling", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  L_a <- (parameters$discrimination != 0) * 1
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, standardize = FALSE, execute = FALSE,
    constraints = list(discrimination = L_a, resid_sd = parameters$resid_params$sd))

  expect_equal(length(samples), 0)
  info <- attr(samples, "spifa_args")
  expect_equal(info$model_type, "cifa")
  expect_equal(info$nobs, nrow(ipixuna_flat))
})

test_that("spifa() surfaces clear errors for malformed constraints/priors", {
  data(ipixuna, package = "spifa")
  parameters <- attr(ipixuna, "parameters")
  nfactors <- ncol(parameters$discrimination)
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  # mismatched constraints$discrimination dimensions
  expect_error(
    spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
          niter = 2, thin = 1, execute = FALSE,
          constraints = list(discrimination = matrix(1, 2, 2))),
    "must be of dimension")

  # mismatched constraints$mgp dimensions (sf data, spatial model)
  expect_error(
    spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
          niter = 2, thin = 1, execute = FALSE,
          constraints = list(mgp = matrix(1, 2, 2))),
    "must be of dimension")

  # malformed priors$easiness$mean length
  expect_error(
    spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
          niter = 2, thin = 1, execute = FALSE,
          priors = list(easiness = list(mean = c(1, 2)))),
    "must be of length")

  # left-hand side of formula must be a matrix-valued column
  expect_error(
    spifa(wealth ~ 1, data = ipixuna_flat, nfactors = nfactors, niter = 2, thin = 1),
    "must be a matrix")
})
