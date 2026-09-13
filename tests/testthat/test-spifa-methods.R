# ---------------------------------------------------------------------------
# summary/print: pure R reshaping of an already-fitted draws array, tested
# once (see test-predict-dic.R for predict()/dic(), which are tested
# against every model_type since they consume freshly-computed,
# model-type-branching C++ output)
# ---------------------------------------------------------------------------

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
  c_block <- as_draws_matrix(subset_draws(samples, variable = "c"))
  expect_equal(smry_burnin$mean, unname(colMeans(c_block[6:10, , drop = FALSE])))

  # a select naming a block absent from this fit (ngp = 0, so no "T") is an
  # error, not a silent zero-row result
  expect_error(summary(samples, select = "T"), "missing")
})

test_that("summary(): burnin >= niter gives a clear error, not a raw seq() crash", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(3, 9), 1] <- 0
  A[c(1, 5, 8), 2] <- 0
  A[c(2, 6, 7, 10), 3] <- 0
  ipixuna_flat <- sf::st_set_geometry(ipixuna, NULL)

  samples <- spifa(items ~ 1, data = ipixuna_flat, nfactors = nfactors,
    niter = 5, thin = 1, constraints = list(discrimination = A))

  # regression test: burnin >= niter used to crash with a raw, confusing
  # "wrong sign in 'by' argument" error from seq(burnin+1, niter, thin);
  # now it hits posterior::subset_draws()'s own clear validation instead
  expect_error(summary(samples, burnin = 5), "iterations")
})

test_that("print.spifa() shows formula, dimensions, and a posterior summary table", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0

  # cifa (a restricted constraints$discrimination): unlike eifa, its Corr is
  # actually estimated, so it has something to show in "Factor model
  # parameters:" below
  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1, constraints = list(discrimination = A))

  out <- capture.output(print(samples))
  out1 <- paste(out, collapse = "\n")
  expect_match(out1, "Formula: items ~ 1", fixed = TRUE)
  expect_match(out1, "100 respondents, 10 items, 3 latent factors", fixed = TRUE)
  expect_match(out1, "c\\[1\\]")
  # Theta/Z/Chol are excluded from the default table (large and/or
  # not directly interpretable)
  expect_no_match(out1, "Theta\\[")
  expect_no_match(out1, "Z\\[")
  expect_no_match(out1, "Chol\\[")

  # summary table is split into item model parameters (c, A) and factor
  # model parameters (B, T, phi, Corr), mirroring brms' grouped summaries
  expect_match(out1, "Item model parameters:", fixed = TRUE)
  expect_match(out1, "Factor model parameters:", fixed = TRUE)
  expect_match(out1, "Corr\\[2,1\\]")
  expect_true(which(out == "Factor model parameters:") >
    which(grepl("^A\\[", out))[1])

  # eifa has no B/T/phi/Corr at all (see the "eifa" test in test-spifa.R):
  # "Factor model parameters:" is omitted entirely rather than printed empty
  samples_eifa <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, thin = 1)
  out_eifa <- paste(capture.output(print(samples_eifa)), collapse = "\n")
  expect_match(out_eifa, "Item model parameters:", fixed = TRUE)
  expect_no_match(out_eifa, "Factor model parameters:")

  # execute = FALSE: no posterior samples to summarise, shouldn't error
  samples_unexecuted <- spifa(items ~ 1, data = ipixuna, nfactors = 3,
    niter = 5, execute = FALSE)
  expect_output(print(samples_unexecuted), "not executed")
})
