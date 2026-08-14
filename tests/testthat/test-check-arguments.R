test_that("check_param_vec fills in defaults and validates length", {
  expect_equal(check_param_vec(list(), "prior_mean", 3, 0), c(0, 0, 0))
  expect_equal(check_param_vec(list(), "prior_mean", 3, c(1, 2, 3)), c(1, 2, 3))
  expect_equal(check_param_vec(list(prior_mean = 5), "prior_mean", 3, 0), c(5, 5, 5))
  expect_equal(check_param_vec(list(prior_mean = c(1, 2, 3)), "prior_mean", 3, 0), c(1, 2, 3))
  expect_error(check_param_vec(list(prior_mean = c(1, 2)), "prior_mean", 3, 0),
               "must be of length 1 or 3")
})

test_that("check_param_mat requires exact-dimension matrices", {
  default <- matrix(0, 2, 3)
  expect_equal(check_param_mat(list(), "A", c(2, 3), default), default)
  m <- matrix(1:6, 2, 3)
  expect_equal(check_param_mat(list(A = m), "A", c(2, 3), default), m)
  expect_error(check_param_mat(list(A = matrix(1, 3, 2)), "A", c(2, 3), default),
               "must be of dimension")
})

test_that("check_param_mat2 accepts a scalar or a matching matrix", {
  expect_equal(check_param_mat2(list(), "B", c(2, 2), 1), matrix(1, 2, 2))
  expect_equal(check_param_mat2(list(B = 5), "B", c(2, 2), 1), matrix(5, 2, 2))
  m <- matrix(1:4, 2, 2)
  expect_equal(check_param_mat2(list(B = m), "B", c(2, 2), 1), m)
  expect_error(check_param_mat2(list(B = matrix(1, 3, 3)), "B", c(2, 2), 1),
               "must be of length 1 or dimension")
})

test_that("check_param_matdiag accepts scalar, vector, or square matrix", {
  expect_equal(check_param_matdiag(list(), "Sigma", 2, diag(2)), diag(2))
  expect_equal(check_param_matdiag(list(Sigma = 3), "Sigma", 2, diag(2)), diag(3, 2, 2))
  expect_equal(check_param_matdiag(list(Sigma = c(1, 2)), "Sigma", 2, diag(2)), diag(c(1, 2)))
  m <- matrix(c(2, 0, 0, 2), 2, 2)
  expect_equal(check_param_matdiag(list(Sigma = m), "Sigma", 2, diag(2)), m)
  expect_error(check_param_matdiag(list(Sigma = c(1, 2, 3)), "Sigma", 2, diag(2)),
               "must be of length")
})
