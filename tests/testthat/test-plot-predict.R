# ---------------------------------------------------------------------------
# plot_predict()
# ---------------------------------------------------------------------------

fit_spatial_samples <- function () {
  data(ipixuna, package = "spifa")
  nitems <- ncol(ipixuna$items)
  nfactors <- 3
  A <- matrix(1, nitems, nfactors)
  A[c(4, 8), 1] <- 0
  A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
  A[c(5, 6), 3] <- 0
  spifa(
    items ~ 1, data = ipixuna, nfactors = nfactors, niter = 5,
    constraints = list(discrimination = A)
  )
}

test_that("plot_predict(): defaults to mean, uses newdata attribute when grid missing", {
  samples <- fit_spatial_samples()
  pred <- predict(samples)

  gg <- plot_predict(pred)
  expect_s3_class(gg, "ggplot")
  expect_equal(nrow(gg$data), nrow(attr(pred, "newdata")) * 3)

  vars <- posterior::variables(pred)
  idx <- grepl(",1\\]$", vars)
  expected <- vapply(seq_len(sum(idx)), function (i) {
    mean(posterior::extract_variable(pred, vars[idx][i]))
  }, numeric(1))
  expect_equal(plot_predict(pred, select = 1)$data$value, expected)
})

test_that("plot_predict(): draws points for points, polygons for polygon newdata", {
  samples <- fit_spatial_samples()
  data(ipixuna, package = "spifa")

  newcoords <- sf::st_make_grid(ipixuna, n = c(5, 1), what = "centers")
  pred_points <- predict(samples, newdata = newcoords)
  gg_points <- plot_predict(pred_points)
  expect_true(inherits(gg_points$layers[[1]]$geom, "GeomSf"))
  expect_true("colour" %in% names(gg_points$layers[[1]]$mapping))

  grid <- sf::st_sf(geometry = sf::st_make_grid(ipixuna, n = c(3, 2)))
  pred_grid <- predict(samples, newdata = grid)
  gg_grid <- plot_predict(pred_grid)
  expect_true("fill" %in% names(gg_grid$layers[[1]]$mapping))
})

test_that("plot_predict(): stat accepts any function of the draws", {
  samples <- fit_spatial_samples()
  data(ipixuna, package = "spifa")
  grid <- sf::st_sf(geometry = sf::st_make_grid(ipixuna, n = c(3, 2)))
  pred <- predict(samples, newdata = grid)
  vars <- posterior::variables(pred)
  idx <- grepl(",1\\]$", vars)

  gg_sd <- plot_predict(pred, stat = sd, select = 1)
  expected_sd <- vapply(seq_len(sum(idx)), function (i) {
    stats::sd(posterior::extract_variable(pred, vars[idx][i]))
  }, numeric(1))
  expect_equal(gg_sd$data$value, expected_sd)

  gg_exc <- plot_predict(pred, stat = function (v) mean(v > 0), select = 1)
  expected_exc <- vapply(seq_len(sum(idx)), function (i) {
    mean(posterior::extract_variable(pred, vars[idx][i]) > 0)
  }, numeric(1))
  expect_equal(gg_exc$data$value, expected_exc)
})

test_that("plot_predict(): rejects a non-function stat", {
  samples <- fit_spatial_samples()
  pred <- predict(samples)
  expect_error(plot_predict(pred, stat = "mean"), "must be a function")
})

test_that("plot_predict(): select restricts which factors are plotted", {
  samples <- fit_spatial_samples()
  data(ipixuna, package = "spifa")
  grid <- sf::st_sf(geometry = sf::st_make_grid(ipixuna, n = c(3, 2)))
  pred <- predict(samples, newdata = grid)

  gg <- plot_predict(pred, select = c(2, 3))
  expect_equal(levels(gg$data$factor), c("theta[2]", "theta[3]"))
  expect_equal(nrow(gg$data), nrow(grid) * 2)
})

test_that("plot_predict(): errors without grid or a newdata attribute", {
  samples <- fit_spatial_samples()
  pred <- predict(samples)
  attr(pred, "newdata") <- NULL

  expect_error(plot_predict(pred), "newdata")
  expect_error(plot_predict(pred, grid = 1:3), "sf/sfc")
})

test_that("plot_predict(): boundary is an extra layer, doesn't affect the data", {
  samples <- fit_spatial_samples()
  data(ipixuna, package = "spifa")
  grid <- sf::st_sf(geometry = sf::st_make_grid(ipixuna, n = c(3, 2)))
  pred <- predict(samples, newdata = grid)
  boundary <- sf::st_convex_hull(sf::st_union(sf::st_geometry(ipixuna)))

  gg_plain <- plot_predict(pred)
  gg_boundary <- plot_predict(pred, boundary = boundary)

  expect_equal(gg_plain$data, gg_boundary$data)
  expect_equal(length(gg_boundary$layers), length(gg_plain$layers) + 1)
})
