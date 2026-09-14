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
    c("c", "A", "Theta", "Z"))

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

test_that("spifa(): a missing item response is kept, a missing predictor value drops the respondent", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 2

  # a single NA item response does not drop the respondent -- the sampler
  # handles it natively as an unobserved auxiliary variable
  ipixuna_item_na <- ipixuna
  items_na <- unclass(ipixuna_item_na$items)
  items_na[1, 1] <- NA
  ipixuna_item_na$items <- items_na

  fit_item_na <- spifa(items ~ 1, data = ipixuna_item_na, nfactors = nfactors,
    ngp = 0, niter = 2, thin = 1)
  expect_equal(attr(fit_item_na, "fit_args")$nobs, nrow(ipixuna))
  expect_true(anyNA(attr(fit_item_na, "fit_args")$response))

  # a single NA predictor value drops that whole respondent (no per-item
  # handling on the predictor side), including from the item responses
  ipixuna_pred_na <- ipixuna
  ipixuna_pred_na$wealth[1] <- NA

  fit_pred_na <- spifa(items ~ wealth, data = ipixuna_pred_na, nfactors = nfactors,
    ngp = 0, niter = 2, thin = 1)
  expect_equal(attr(fit_pred_na, "fit_args")$nobs, nrow(ipixuna) - 1)

  # both at once, on different respondents: only the predictor-NA
  # respondent is dropped, the item-NA respondent is kept (with its NA)
  ipixuna_both_na <- ipixuna_pred_na
  items_na2 <- unclass(ipixuna_both_na$items)
  items_na2[2, 1] <- NA
  ipixuna_both_na$items <- items_na2

  fit_both_na <- spifa(items ~ wealth, data = ipixuna_both_na, nfactors = nfactors,
    ngp = 0, niter = 2, thin = 1)
  expect_equal(attr(fit_both_na, "fit_args")$nobs, nrow(ipixuna) - 1)
  expect_true(anyNA(attr(fit_both_na, "fit_args")$response))

  # transformed predictor terms (e.g. poly()) are still matched correctly
  # when checking predictor completeness
  fit_poly <- spifa(items ~ poly(wealth, 2), data = ipixuna, nfactors = nfactors,
    ngp = 0, niter = 2, thin = 1)
  expect_equal(attr(fit_poly, "fit_args")$nobs, nrow(ipixuna))
})

# ---------------------------------------------------------------------------
# update.spifa()
# ---------------------------------------------------------------------------

test_that("update.spifa() warm-starts from the last draw, across model types", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0

  # warm-start values must match the previous fit's last draw exactly (a
  # single subsequent Gibbs step would already move them, so check the
  # *initial* values actually passed to spifa_cpp(), not the first new draw)
  expect_warm_started <- function (samples, blocks = c("c", "A")) {
    samples2 <- update(samples, niter = 3)
    for (block in blocks) {
      last <- as.numeric(as_draws_matrix(subset_draws(samples, variable = block))[
        niterations(samples), ])
      used <- if (block == "c") {
        attr(samples2, "fit_args")$c_initial
      } else if (block == "A") {
        as.numeric(attr(samples2, "fit_args")$A_initial)
      }
      expect_equal(unname(used), last)
    }
    expect_equal(niterations(samples2), 3)
    samples2
  }

  # eifa: Corr/Chol absent -- reused fixed value, not "warm-started"
  eifa <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 10, thin = 1)
  eifa2 <- expect_warm_started(eifa)
  expect_equal(attr(eifa2, "fit_args")$R_initial, attr(eifa, "fit_args")$R_initial)
  expect_setequal(variables(eifa2, with_indices = FALSE), c("c", "A", "Theta", "Z"))

  # cifa: Corr/Chol present and genuinely warm-started
  cifa <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 10, thin = 1, constraints = list(discrimination = A))
  cifa2 <- expect_warm_started(cifa)
  last_corr <- unname(as_draws_matrix(subset_draws(cifa, variable = "Corr"))[10, ])
  expect_true(any(diag(attr(cifa2, "fit_args")$R_initial) == 1))
  expect_false(isTRUE(all.equal(attr(cifa2, "fit_args")$R_initial, diag(nfactors))))

  # cifa_pred: B present
  cifa_pred <- spifa(items ~ wealth, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 10, thin = 1, constraints = list(discrimination = A))
  cifa_pred2 <- expect_warm_started(cifa_pred)
  expect_equal(dim(attr(cifa_pred2, "fit_args")$B_initial), c(1, nfactors))

  # spifa: T/phi present
  spifa_fit <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, niter = 10,
    thin = 1, constraints = list(discrimination = A, loading = diag(nfactors)),
    priors = list(
      loading = list(initial = 0.6, mean = 0.6, sd = 0.4),
      range = list(initial = 200, mean = 200, sd = 0.4)))
  spifa2 <- expect_warm_started(spifa_fit)
  expect_length(attr(spifa2, "fit_args")$sigmas_gp_initial, nfactors)
  expect_length(attr(spifa2, "fit_args")$phi_gp_initial, nfactors)
})

test_that("update.spifa() errors clearly for an unexecuted fit", {
  data(ipixuna, package = "spifa")
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, niter = 5,
    execute = FALSE)
  expect_error(update(samples, niter = 5), "not executed")
})

test_that("update.spifa() resumes the adaptive-MH proposal instead of restarting it", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0

  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 50, thin = 1, constraints = list(discrimination = A))
  mcmc_state <- attr(samples, "mcmc_state")
  expect_equal(dim(mcmc_state$adap_Sigma), c(3, 3))
  expect_type(mcmc_state$adap_scale, "double")

  samples2 <- update(samples, niter = 10)
  # the continuation is warm-started from the previous run's ending
  # adaptive state, not object's own *original* adap_Sigma/adap_scale
  expect_equal(attr(samples2, "fit_args")$adap_Sigma, mcmc_state$adap_Sigma)
  expect_equal(attr(samples2, "fit_args")$adap_scale, mcmc_state$adap_scale)

  # and samples2 carries its own (further-adapted) ending state forward, so
  # a chain of update() calls keeps refining rather than resetting each time
  mcmc_state2 <- attr(samples2, "mcmc_state")
  expect_false(isTRUE(all.equal(mcmc_state2$adap_Sigma, mcmc_state$adap_Sigma)))
})
