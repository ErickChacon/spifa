# ---------------------------------------------------------------------------
# summary/print: pure R reshaping of an already-fitted draws array
# ---------------------------------------------------------------------------

test_that("summary.spifa(): returns posterior summaries with burnin/thin/select", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(1, 3, 5), 1] <- 0
  A[c(2, 6, 9), 2] <- 0
  A[c(4, 7, 10), 3] <- 0

  samples <- spifa(
    items ~ 1, data = ipixun, nfactors = nfactors, ngp = 0, niter = 10,
    constraints = list(discrimination = A)
  )

  smry <- summary(samples, select = "c")
  expect_true(all(c("variable", "q2.5", "median", "q97.5") %in% names(smry)))
  expect_equal(nrow(smry), ncol(ipixun$items))

  smry_burnin <- summary(samples, burnin = 5, select = "c")
  c_block <- as_draws_matrix(subset_draws(samples, variable = "c"))
  expect_equal(smry_burnin$mean, unname(colMeans(c_block[6:10, , drop = FALSE])))

  # error when selecting an absent term
  expect_error(summary(samples, select = "T"), "missing")
})

test_that("summary.spifa()/print.spifa(): exclude restricted A parameters", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(1, 3, 5), 1] <- 0
  A[c(2, 6, 9), 2] <- 0
  A[c(4, 7, 10), 3] <- 0

  samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, constraints = list(discrimination = A))

  smry <- summary(samples, select = "A")
  expect_false("A[1,1]" %in% smry$variable)
  expect_true("A[2,1]" %in% smry$variable)
  expect_equal(nrow(smry), sum(A == 1))

  out <- paste(capture.output(print(samples)), collapse = "\n")
  expect_no_match(out, "A\\[1,1\\]")
  expect_match(out, "A\\[2,1\\]")

  # eifa's default lower-triangular restriction is excluded the same way
  samples_eifa <- spifa(
    items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0, niter = 5
  )
  smry_eifa <- summary(samples_eifa, select = "A")
  expect_equal(nrow(smry_eifa),
    sum(lower.tri(matrix(NA, nitems, nfactors), diag = TRUE)))
})

test_that("summary(): burnin >= niter gives an error", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(3, 9), 1] <- 0
  A[c(1, 5, 8), 2] <- 0
  A[c(2, 6, 7, 10), 3] <- 0

  samples <- spifa(
    items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0, niter = 5,
    constraints = list(discrimination = A)
  )

  expect_error(summary(samples, burnin = 5), "iterations")
})

test_that("print.spifa(): shows formula, dimensions, and summary table", {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0

  # cifa
  samples <- spifa(
    items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0,
    niter = 5, constraints = list(discrimination = A)
  )

  out <- capture.output(print(samples))
  out1 <- paste(out, collapse = "\n")
  expect_match(out1, "Formula: items ~ 1", fixed = TRUE)
  expect_match(out1, "100 respondents, 10 items, 3 latent factors", fixed = TRUE)
  expect_match(out1, "c\\[1\\]")
  # Theta/Z/Chol are excluded from the default table
  expect_no_match(out1, "Theta\\[")
  expect_no_match(out1, "Z\\[")
  expect_no_match(out1, "Chol\\[")

  # summary splitted into item (c,A) and factor (B, T, phi, Corr) parameters
  expect_match(out1, "Item model parameters:", fixed = TRUE)
  expect_match(out1, "Factor model parameters:", fixed = TRUE)
  expect_match(out1, "Corr\\[2,1\\]")
  expect_true(which(out == "Factor model parameters:") >
    which(grepl("^A\\[", out))[1])

  # eifa: has no B/T/phi/Corr
  samples_eifa <- spifa(
    items ~ 1, data = ipixuna, nfactors = nfactors, ngp = 0, niter = 5
  )
  # "Factor model parameters:" is omitted entirely rather than printed empty
  out_eifa <- paste(capture.output(print(samples_eifa)), collapse = "\n")
  expect_match(out_eifa, "Item model parameters:", fixed = TRUE)
  expect_no_match(out_eifa, "Factor model parameters:")

  # execute = FALSE: no posterior samples to summarise, shouldn't error
  samples_unexecuted <- spifa(
    items ~ 1, data = ipixuna, nfactors = 3, niter = 5, execute = FALSE
  )
  expect_output(print(samples_unexecuted), "not executed")
})
