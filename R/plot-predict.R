# clean map look for plot_predict(): no axis titles, a light dashed
# reference grid drawn on top of the fill, small legend, no strip background
theme_map <- function (base_size = 9, legend = "bottom") {
  theme_bw(base_size = base_size) +
    theme(axis.title = element_blank(),
          panel.grid = element_line(linetype = 2, colour = "grey80", linewidth = 0.25),
          panel.background = element_blank(),
          panel.ontop = TRUE,
          legend.position = legend,
          legend.title = element_blank(),
          legend.key.height = unit(0.4, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(size = rel(1)))
}

#' @title Map Posterior Predictions from a Fitted spifa Model
#'
#' @description
#' Maps a posterior summary of each latent factor (\code{stat}: the
#' posterior mean by default, or any other function of the draws) at the
#' locations \code{pred} was predicted at. Draws points if those locations
#' are points (e.g. the observed households), or a filled polygon map if
#' they are polygons (e.g. a prediction grid).
#'
#' @details
#' \code{pred} is computed once via \code{\link{predict.spifa}} and can be
#' reused across multiple \code{plot_predict()} calls -- e.g. different
#' \code{stat}/\code{select} values -- without calling \code{predict()}
#' again, which is normally the expensive step.
#'
#' @param pred A \code{\link[posterior]{draws_array}} as returned by
#' \code{\link{predict.spifa}} for a spatial model.
#' @param grid An \code{sf}/\code{sfc} object of the locations \code{pred}
#' was predicted at (points or polygons), in the same row order as
#' \code{pred}'s locations. Defaults to \code{attr(pred, "newdata")},
#' attached automatically by \code{\link{predict.spifa}}; only needed
#' explicitly if that attribute is missing.
#' @param select Factors to plot: an integer vector of factor indices, or
#' \code{NULL} (default) for all.
#' @param stat A function of one numeric vector (the draws at one
#' location/factor) returning a single number -- e.g. \code{mean} (the
#' default), \code{sd}, \code{median}, or \code{function(v) mean(v > 1)}
#' for an exceedance probability.
#' @param boundary Optional \code{sf}/\code{sfc} polygon(s) drawn as an
#' outline (e.g. the prediction area) on top of the map; has no effect on
#' the computed/plotted values.
#' @param facet_scales The \code{scales} argument of
#' \code{\link[ggplot2]{facet_wrap}}.
#' @param ncol Number of columns in the facet grid; default lets
#' \code{\link[ggplot2]{facet_wrap}} choose.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object with a clean map theme already applied
#' (no axis titles, light dashed reference grid, small bottom legend). A
#' fill (polygons) or colour (points) scale can be added as usual; a
#' \code{\link[ggplot2]{theme}}/\code{theme_*()} added afterwards overrides
#' it as normal.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' library(sf)
#' data(ipixuna)
#'
#' nitems <- ncol(ipixuna$items)
#' nfactors <- 3
#' A <- matrix(1, nitems, nfactors)
#' A[c(4, 8), 1] <- 0
#' A[c(2, 4, 5, 6, 7, 8, 10), 2] <- 0
#' A[c(5, 6), 3] <- 0
#'
#' samples <- spifa(items ~ 1, data = ipixuna, nfactors = nfactors, niter = 5,
#'   constraints = list(discrimination = A))
#'
#' # at the observed household locations (points)
#' pred <- predict(samples)
#' plot_predict(pred)
#'
#' # on a grid (polygons); predict() takes the cell centroids for the
#' # spatial kernel but keeps the polygons for plot_predict() to map
#' grid <- st_sf(geometry = st_make_grid(ipixuna, n = c(3, 2)))
#' pred_grid <- predict(samples, newdata = grid)
#' plot_predict(pred_grid)
#' plot_predict(pred_grid, stat = sd)
#' plot_predict(pred_grid, stat = function (v) mean(v > 0)) # exceedance
#' }
#'
#' @export
plot_predict <- function (pred, grid = NULL, select = NULL, stat = mean,
                           boundary = NULL, facet_scales = "fixed", ncol = NULL, ...) {
  if (!is.function(stat)) {
    stop("`stat` must be a function (e.g. `mean`, `sd`, or ",
         "`function(v) mean(v > 1)`).", call. = FALSE)
  }
  if (is.null(grid)) {
    grid <- attr(pred, "newdata")
  }
  if (is.null(grid)) {
    stop("`grid` was not supplied and `pred` has no \"newdata\" attribute ",
         "(only attached by predict.spifa() for spatial models) -- pass it explicitly.",
         call. = FALSE)
  }
  if (!(inherits(grid, "sf") || inherits(grid, "sfc"))) {
    stop("`grid` must be an sf/sfc object.", call. = FALSE)
  }

  vars <- posterior::variables(pred)
  factor_idx <- as.integer(sub(".*,(\\d+)\\]$", "\\1", vars))
  if (is.null(select)) {
    select <- sort(unique(factor_idx))
  }

  summ <- posterior::summarise_draws(pred, value = stat)

  geom <- sf::st_geometry(grid)
  df <- do.call(rbind, lapply(select, function (j) {
    data.frame(factor = paste0("theta[", j, "]"), value = summ$value[factor_idx == j])
  }))
  df$factor <- factor(df$factor, levels = paste0("theta[", select, "]"))
  df <- sf::st_sf(df, geometry = rep(geom, length(select)))

  geom_type <- unique(sf::st_geometry_type(geom))
  is_point <- all(geom_type %in% c("POINT", "MULTIPOINT"))

  gg <- ggplot(df)
  if (is_point) {
    gg <- gg + geom_sf(aes(colour = value))
  } else {
    gg <- gg + geom_sf(aes(fill = value), colour = NA)
  }
  if (!is.null(boundary)) {
    gg <- gg + geom_sf(data = boundary, fill = NA, linewidth = 0.5, colour = "grey50")
  }
  gg +
    facet_wrap(~ factor, scales = facet_scales, ncol = ncol, labeller = label_parsed) +
    theme_map()
}
