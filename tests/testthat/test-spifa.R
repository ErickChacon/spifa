# ---------------------------------------------------------------------------
# Per-model-type: eifa, cifa, cifa_pred, spifa, spifa_pred
# ---------------------------------------------------------------------------

test_that("eifa", {
  data(ipixuna, package = "spifa")
  nfactors <- 3
  nitems <- ncol(ipixuna$items)

  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1)

  expect_equal(attr(samples, "fit_args")$model_type, "eifa")
  expect_equal(attr(samples, "fit_args")$constrain_L,
    lower.tri(matrix(NA, nitems, nfactors), diag = TRUE) * 1)
  expect_equal(length(attr(samples, "fit_args")$response), nrow(ipixuna) * nitems)
  expect_setequal(variables(samples, with_indices = FALSE),
    c("c", "A", "Theta", "Z", "Chol", "Corr"))

  # standardize does not affect eifa models
  set.seed(42)
  samples_raw <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
     niter = 5, thin = 1, standardize = FALSE)
  set.seed(42)
  samples_std <- spifa(items ~ 1, data = sf::st_set_geometry(ipixuna, NULL),
    nfactors = nfactors, niter = 5, thin = 1, standardize = TRUE)
  expect_equal(as_draws_matrix(samples_raw),
    as_draws_matrix(samples_std), ignore_attr = "fit_args")
})

test_that("cifa", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(1, 3, 5), 1] <- 0
  A[c(2, 6, 9), 2] <- 0
  A[c(4, 7, 10), 3] <- 0

  samples <- spifa(
    items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1,
    constraints = list(discrimination = A))

  expect_equal(attr(samples, "fit_args")$model_type, "cifa")
  expect_setequal(variables(samples, with_indices = FALSE),
    c("c", "A", "Theta", "Z", "Chol", "Corr"))

  # verify zero loadings based on constraints
  a_samples <- as_draws_matrix(subset_draws(samples, variable = "A"))
  idx <- as.integer(sub("A\\[(\\d+),(\\d+)\\]", "\\1", colnames(a_samples)))
  jdx <- as.integer(sub("A\\[(\\d+),(\\d+)\\]", "\\2", colnames(a_samples)))
  zero_mask <- A[cbind(idx, jdx)] == 0
  expect_true(all(a_samples[, zero_mask] == 0))
})

test_that("cifa_pred", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(2, 7), 1] <- 0
  A[c(1, 4, 9, 10), 2] <- 0
  A[c(3, 6), 3] <- 0

  samples <- spifa(
    items ~ poly(wealth, 2), data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1,
    constraints = list(discrimination = A))

  expect_equal(attr(samples, "fit_args")$model_type, "cifa_pred")
  expect_equal(dim(attr(samples, "fit_args")$predictors), c(nrow(ipixuna), 2))
  expect_setequal(variables(samples, with_indices = FALSE),
    c("c", "A", "Theta", "Z", "Chol", "Corr", "B"))
  b_block <- subset_draws(samples, variable = "B")
  expect_equal(ndraws(b_block), 5)
  expect_equal(nvariables(b_block), 2 * nfactors)

  # standardize rescales theta (latent abilities) to unit variance
  set.seed(1)
  samples_std <- spifa(
    items ~ wealth, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 200, thin = 1, standardize = TRUE,
    constraints = list(discrimination = A))
  theta <- as_draws_matrix(subset_draws(samples_std, variable = "Theta"))
  factor_idx <- as.integer(sub(".*,(\\d+)\\]$", "\\1", colnames(theta)))
  for (k in seq_len(nfactors)) {
    expect_true(abs(sd(as.numeric(theta[, factor_idx == k])) - 1) < 0.1)
  }
})

test_that("spifa", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0

  # default case
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
    niter = 5, thin = 1,
    constraints = list(discrimination = A))

  expect_equal(attr(samples, "fit_args")$model_type, "spifa")
  expect_equal(attr(samples, "fit_args")$ngp, nfactors)
  expect_equal(attr(samples, "fit_args")$constrain_T, diag(1, nfactors, nfactors))
  expect_equal(attr(samples, "predict_setup")$coordinates, sf::st_geometry(ipixuna))
  expect_setequal(variables(samples, with_indices = FALSE),
    c("c", "A", "Theta", "Z", "Chol", "Corr", "T", "phi"))
  expect_equal(nvariables(subset_draws(samples, variable = "T")), nfactors)
  expect_equal(nvariables(subset_draws(samples, variable = "phi")), nfactors)

  # custom gps and factors relationship
  fit_spifa <- function (loading, ngp) {
    spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = ngp,
      niter = 5, thin = 1,
      constraints = list(discrimination = A, loading = loading),
      priors = list(
        loading = list(initial = rep(0.6, sum(loading)), mean = rep(0.6, sum(loading)), sd = 0.4),
        range = list(initial = rep(200, ngp), mean = rep(200, ngp), sd = 0.3))
    )
  }

  # 1 GP
  loading_shared <- matrix(c(1, 1, 0), nfactors, 1)
  samples_shared <- fit_spifa(loading_shared, 1)
  expect_equal(nvariables(subset_draws(samples_shared, variable = "T")),
    sum(loading_shared))
  expect_equal(nvariables(subset_draws(samples_shared, variable = "phi")), 1)

  # 2 GP
  loading_partial <- matrix(c(1, 0, 0, 0, 1, 1), nfactors, 2)
  samples_partial <- fit_spifa(loading_partial, 2)
  expect_equal(nvariables(subset_draws(samples_partial, variable = "T")),
    sum(loading_partial))
  expect_equal(nvariables(subset_draws(samples_partial, variable = "phi")), 2)
})

test_that("spifa(): standardize = TRUE (the default) works with no predictors", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0

  # regression test: standardize = TRUE used to crash for a plain spatial
  # model with no predictors ("Mat::submat(): indices or size out of
  # bounds") -- betas_samples has 0 rows when p (predictor count) is 0,
  # but the C++ standardize step indexed into it unconditionally
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
    niter = 20, thin = 1,
    constraints = list(discrimination = A))

  expect_equal(attr(samples, "fit_args")$model_type, "spifa")
})

test_that("spifa_pred", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0

  samples <- spifa(
    items ~ wealth, data = ipixuna, nfactors = nfactors, niter = 5, thin = 1,
    constraints = list(discrimination = A, loading = diag(nfactors)),
    priors = list(
      loading = list(initial = 0.6, mean = 0.6, sd = 0.4),
      range = list(initial = 200, mean = 200, sd = 0.4)))

  expect_equal(attr(samples, "fit_args")$model_type, "spifa_pred")
  expect_equal(dim(attr(samples, "fit_args")$predictors), c(nrow(ipixuna), 1))
  expect_setequal(variables(samples, with_indices = FALSE),
    c("c", "A", "Theta", "Z", "Chol", "Corr", "T", "phi", "B"))
})

# ---------------------------------------------------------------------------
# Sampler mechanics and argument validation
# ---------------------------------------------------------------------------

test_that("spifa(): burnin, thin, niter", {
  data(ipixuna, package = "spifa")
  nfactors <- 2

  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1, burnin = 0)
  expect_equal(dim(samples)[1], 5)

  samples_burnin <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1, burnin = 10)
  expect_equal(dim(samples_burnin)[1], 5)

  # non-divisible niter and thin
  samples_edge <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 8, thin = 3, burnin = 2)
  expect_equal(dim(samples_edge)[1], ceiling(8 / 3))
  expect_equal(attr(samples_edge, "fit_args")$niter, 7)
})

test_that("spifa(): execute = FALSE", {
  data(ipixuna, package = "spifa")
  nfactors <- 2
  L <- matrix(1, ncol(ipixuna$items), nfactors)
  L[1:3, nfactors] <- 0
  easiness_mean <- rep(0.5, ncol(ipixuna$items))

  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
    niter = 5, thin = 1, execute = FALSE,
    constraints = list(discrimination = L),
    priors = list(easiness = list(mean = easiness_mean)))

  expect_equal(length(samples), 0)
  info <- attr(samples, "fit_args")
  expect_equal(info$model_type, "spifa")
  expect_equal(info$nobs, nrow(ipixuna))
  expect_equal(info$constrain_L, L)
  expect_equal(info$c_prior_mean, easiness_mean)
  expect_equal(attr(samples, "predict_setup")$coordinates, sf::st_geometry(ipixuna))
})

test_that("spifa(): malformed constraints/priors", {
  data(ipixuna, package = "spifa")
  nfactors <- 2

  # mismatched constraints$discrimination dimensions
  expect_error(
    spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
          niter = 2, thin = 1, execute = FALSE,
          constraints = list(discrimination = matrix(1, 2, 2))),
    "must be of dimension")

  # mismatched constraints$loading dimensions (spatial model)
  expect_error(
    spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
          niter = 2, thin = 1, execute = FALSE,
          constraints = list(loading = matrix(1, 3, 3))),
    "must be of dimension")

  # malformed priors$easiness$mean length
  expect_error(
    spifa(items ~ 1, data = ipixuna, nfactors = nfactors,
          niter = 2, thin = 1, execute = FALSE,
          priors = list(easiness = list(mean = c(1, 2)))),
    "must be of length")

  # left-hand side of formula must be a matrix-valued column
  expect_error(
    spifa(wealth ~ 1, data = ipixuna, nfactors = nfactors, niter = 2, thin = 1),
    "must be a matrix")
})
