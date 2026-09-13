# remove restricted discrimination parameters
drop_restricted <- function (x) {

  fit_args <- attr(x, "fit_args")
  predict_setup <- attr(x, "predict_setup")

  index <- which(fit_args$constrain_L == 0, arr.ind = TRUE)
  a_restricted <- paste0("A[", index[, 1], ",", index[, 2], "]")
  a_unrestricted <- setdiff(posterior::variables(x), a_restricted)
  x <- posterior::subset_draws(x, variable = a_unrestricted)

  attr(x, "fit_args") <- fit_args
  attr(x, "predict_setup") <- predict_setup
  return(x)
}

# plotmath parse of parameters
parse_parameter <- function (parameter) {
  to_math <- function (p) gsub("\\[(.+),(.+)\\]", "[list(\\1,\\2)]", p)
  factor(to_math(as.character(parameter)), levels = to_math(levels(parameter)))
}

# converts samples to long format
draws_long <- function (x) {
  wide <- x |> posterior::as_draws_df() |> tibble::as_tibble()
  pars <- setdiff(names(wide), c(".chain", ".iteration", ".draw"))
  wide |>
    dplyr::mutate(iteration = .iteration) |>
    tidyr::pivot_longer(dplyr::all_of(pars), names_to = "parameter", values_to = "value") |>
    dplyr::mutate(parameter = factor(parameter, levels = pars)) |>
    dplyr::select(iteration, parameter, value)
}

theme_trace <- function (legend = "none") {
  theme_minimal() +
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linetype = "dashed"),
          axis.line = element_line(),
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0),
          legend.position = legend)
}

#' @title Traceplot of Samples
#'
#' @description
#' Draws MCMC traceplots (iteration vs. value) directly from a fitted
#' \code{spifa} model, in one of two formats: faceted (\code{facet = TRUE},
#' the default; one panel per parameter, with its own free y-scale, so a
#' slow-mixing or small-variance parameter isn't visually flattened by
#' others sharing the same axis) or overlaid (\code{facet = FALSE}; all
#' series on a single panel, for a quick glance at overall convergence).
#' Reshapes draws directly, so it works on \code{x} with no intermediate
#' conversion.
#'
#' @param x A fitted \code{spifa} model.
#' @param select Parameters to plot, passed to the \code{variable} argument
#' of \code{\link[posterior]{subset_draws}}: either a block name (e.g.
#' \code{"A"}, matching every parameter in that block) or one or more full
#' indexed names (e.g. \code{"c[1]"}, \code{paste0("A[", 1:10, ",1]")}).
#' @param facet Logical; if \code{TRUE} (default), draw one panel per
#' parameter; if \code{FALSE}, overlay every series on a single panel. See
#' Description.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after \code{burnin}.
#' @param nshow If \code{select} matches more than \code{nshow} parameters,
#' a random (sorted) subsample of \code{nshow} of them is shown instead of
#' all of them. Set to \code{NULL} to always show every matched parameter.
#' @param ncol Number of columns in the facet grid (\code{facet = TRUE}
#' only); defaults to a single column.
#' @param legend Legend position (\code{facet = FALSE} only): one of
#' \code{"bottom"}, \code{"right"}, \code{"none"}, or a logical. Defaults to
#' \code{NULL}, which shows the legend unless there are more than 10
#' series, since a legend that large stops being readable and can dwarf the
#' plot itself.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 1000)
#'
#' plot_trace(samples, select = "c", facet = TRUE)
#' plot_trace(samples, select = "c", facet = FALSE)
#'
#' # more than nshow (10) parameters: a random subsample is shown
#' plot_trace(samples, select = "A", facet = TRUE, nshow = 6)
#'
#' # explicit selection instead of a random subsample
#' plot_trace(samples, select = paste0("A[", 1:10, ",1]"), facet = TRUE)
#' }
#'
#' @export
plot_trace <- function (x, select, facet = TRUE,
                        burnin = 0, thin = 1, nshow = 10, ncol = 1,
                        legend = NULL, ...) {
  x <- drop_restricted(x)

  # create data
  niter <- posterior::niterations(x)
  df <- x |>
    posterior::subset_draws(variable = select, iteration = (burnin + 1):niter) |>
    posterior::thin_draws(thin) |>
    draws_long() |>
    dplyr::mutate(parameter = parse_parameter(parameter))

  # random sorted subsample when there are more parameters than nshow
  pars <- unique(df$parameter)
  if (!is.null(nshow) && length(pars) > nshow) {
    pars_keep <- sort(sample(pars, nshow))
    df <- dplyr::filter(df, parameter %in% pars_keep)
  }

  # two figure types
  if (!facet) {
    # remove legend if more than 10 parameters
    if (is.null(legend)) {
      legend <- if (length(unique(df$parameter)) > 10) "none" else "bottom"
    } else if (is.logical(legend)) {
      legend <- if (legend) "bottom" else "none"
    }
    gg <- ggplot(df, aes(iteration, value, group = parameter, col = parameter)) +
      geom_path(alpha = 0.6) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_colour_discrete(labels = function(x) parse(text = x)) +
      labs(x = "Iteration", y = "Value", col = "Parameter") +
      theme_trace(legend = legend)
  } else {
    gg <- ggplot(df, aes(iteration, value, col = parameter)) +
      geom_path(linewidth = 0.2) +
      facet_wrap(~ parameter, ncol = ncol, scales = "free_y", strip.position = "right",
                 labeller = label_parsed) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(x = "Iteration", y = "Value") +
      theme_trace(legend = "none")
  }

  return(gg)
}

#' @title Density Plot of Samples
#'
#' @description
#' Draws posterior density curves directly from a fitted \code{spifa}
#' model, in one of two formats: overlaid (\code{facet = FALSE}, the
#' default; stacked, slightly-overlapping ridgeline plot via
#' \code{\link[ggridges]{geom_density_ridges}} -- the more insightful
#' default for comparing many parameters' shapes and locations at a glance) or
#' faceted (\code{facet = TRUE}; one panel per parameter). Shows the
#' density shape only -- for credible intervals and point estimates, see
#' \code{\link{plot_interval}}. The density itself is computed directly by
#' \code{\link[ggplot2]{geom_density}}/\code{\link[ggridges]{geom_density_ridges}}
#' from the raw draws, reshaped via \code{\link[posterior]{as_draws_df}}.
#'
#' @param x A fitted \code{spifa} model.
#' @param select Parameters to plot, as in \code{\link{plot_trace}}.
#' @param facet Logical; if \code{TRUE}, draw one panel per parameter; if
#' \code{FALSE} (default), overlay every parameter as a ridgeline plot. See
#' Description.
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after \code{burnin}.
#' @param nshow As in \code{\link{plot_trace}}.
#' @param ncol Number of columns in the facet grid (\code{facet = TRUE}
#' only); defaults to a single column.
#' @param facet_scales The \code{scales} argument of
#' \code{\link[ggplot2]{facet_wrap}} (\code{facet = TRUE} only): one of
#' \code{"free"} (default; each panel gets its own x/y-axis), \code{"free_y"}
#' (one shared x-axis, line shown only on the bottom panel of each column,
#' matching \code{\link{plot_trace}}'s facet -- useful for comparing
#' parameters on a similar scale), \code{"free_x"}, or \code{"fixed"}.
#' Parameters with very different value ranges can render compressed/hard
#' to read with a shared x-axis.
#' @param scale Amount of vertical overlap between ridges (\code{facet =
#' FALSE} only), passed to \code{\link[ggridges]{geom_density_ridges}}.
#' Defaults to \code{1.2}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 1000)
#'
#' plot_density(samples, select = "c", facet = FALSE)
#' plot_density(samples, select = "c", facet = TRUE)
#'
#' # more than nshow (10) parameters: a random subsample is shown
#' plot_density(samples, select = "A", facet = FALSE, nshow = 6)
#' }
#'
#' @export
plot_density <- function (x, select, facet = FALSE,
                           burnin = 0, thin = 1, nshow = 10, ncol = 1,
                           facet_scales = "free", scale = 1.2, ...) {
  x <- drop_restricted(x)

  # create data
  niter <- posterior::niterations(x)
  df <- x |>
    posterior::subset_draws(variable = select, iteration = (burnin + 1):niter) |>
    posterior::thin_draws(thin) |>
    draws_long() |>
    dplyr::mutate(parameter = parse_parameter(parameter))

  # random sorted subsample when there are more parameters than nshow
  pars <- unique(df$parameter)
  if (!is.null(nshow) && length(pars) > nshow) {
    pars_keep <- sort(sample(pars, nshow))
    df <- dplyr::filter(df, parameter %in% pars_keep)
  }

  # figure types
  if (!facet) {
    gg <- ggplot(df, aes(value, parameter, fill = parameter)) +
      ggridges::geom_density_ridges(alpha = 0.5, scale = scale, linewidth = 0.4,
                                     rel_min_height = 0.01) +
      scale_y_discrete(labels = function(x) parse(text = x)) +
      labs(x = "Value", y = NULL) +
      theme_trace(legend = "none")
  } else {
    gg <- ggplot(df, aes(value, fill = parameter)) +
      geom_density(alpha = 0.5, linewidth = 0.4, trim = TRUE) +
      facet_wrap(~ parameter, ncol = ncol, scales = facet_scales,
                 strip.position = "right", labeller = label_parsed) +
      labs(x = "Value", y = "Density") +
      theme_trace(legend = "none")
  }

  return(gg)
}

#' @title Interval Plot of Samples
#'
#' @description
#' Draws a caterpillar/forest plot of posterior credible intervals directly
#' from a fitted \code{spifa} model: one row per parameter, with a thin
#' line for the \code{prob_outer} interval, a thick line for the
#' \code{prob} interval, and a point at the \code{point_est}, computed
#' directly from the raw draws (quantiles/median/mean). Unlike
#' \code{\link{plot_trace}}/\code{\link{plot_density}}, this is a single
#' combined view by design -- comparing intervals side by side is the whole
#' point, so there is no faceted alternative -- but \code{sort} can reorder
#' parameters by their point estimate, which a facet can't do usefully
#' across independent panels.
#'
#' @param x A fitted \code{spifa} model.
#' @param select Parameters to plot, as in \code{\link{plot_trace}}.
#' @param horizontal Logical; if \code{FALSE} (default), parameters run
#' along the x-axis and values run along the y-axis, matching the
#' \code{ci_intervals()} convention used in the SPIFA paper's own figures;
#' if \code{TRUE}, the axes are swapped (a forest-plot layout).
#' @param burnin Number of initial iterations to discard.
#' @param thin Thinning interval applied after \code{burnin}.
#' @param nshow As in \code{\link{plot_trace}} (a random subsample when
#' \code{select} matches more than \code{nshow} parameters), but
#' \code{NULL} (show every matched parameter) by default: unlike a faceted
#' plot, a single interval plot stays readable with many more than 10 rows.
#' @param prob Width of the thick (inner) credible interval (a central
#' quantile interval). Defaults to \code{0.5}.
#' @param prob_outer Width of the thin (outer) credible interval. Defaults
#' to \code{0.9}.
#' @param point_est Either \code{"median"} (default) or \code{"mean"}.
#' @param sort Logical; if \code{TRUE}, reorder parameters by their point
#' estimate instead of their natural order. Defaults to \code{FALSE}.
#' @param reference Optional reference values to overlay (e.g. the true
#' values in a simulation study), as a fourth marker alongside the interval
#' and point estimate. Only valid when \code{select} is a single block name
#' (e.g. \code{"A"}): an unnamed vector or matrix matching that block's own
#' shape (e.g. \code{parameters$discrimination}, an \code{nitems x nfactors}
#' matrix). Structurally-restricted parameters (dropped internally before
#' plotting) are silently ignored if present in \code{reference}.
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' \donttest{
#' data(ipixuna)
#' samples <- spifa(items ~ 1, data = ipixuna, nfactors = 3, ngp = 0, niter = 1000)
#'
#' plot_interval(samples, select = "c")
#' plot_interval(samples, select = "c", sort = TRUE)
#' plot_interval(samples, select = "A", horizontal = TRUE)
#'
#' # overlay the true simulated discrimination values
#' parameters <- attr(ipixuna, "parameters")
#' plot_interval(samples, select = "A", reference = parameters$discrimination)
#' }
#'
#' @export
plot_interval <- function (x, select, horizontal = FALSE,
                            burnin = 0, thin = 1, nshow = NULL,
                            prob = 0.5, prob_outer = 0.9,
                            point_est = c("median", "mean"), sort = FALSE,
                            reference = NULL, ...) {
  x <- drop_restricted(x)
  point_est <- match.arg(point_est)

  # create data: one row per parameter
  niter <- posterior::niterations(x)
  df <- x |>
    posterior::subset_draws(variable = select, iteration = (burnin + 1):niter) |>
    posterior::thin_draws(thin) |>
    draws_long() |>
    dplyr::group_by(parameter) |>
    dplyr::summarise(
      ll = stats::quantile(value, (1 - prob_outer) / 2, names = FALSE),
      l = stats::quantile(value, (1 - prob) / 2, names = FALSE),
      m = if (point_est == "median") stats::median(value) else mean(value),
      h = stats::quantile(value, 1 - (1 - prob) / 2, names = FALSE),
      hh = stats::quantile(value, 1 - (1 - prob_outer) / 2, names = FALSE),
      .groups = "drop")

  # add reference data frame
  if (!is.null(reference)) {
    if (length(select) != 1 || grepl("[", select, fixed = TRUE)) {
      stop("`reference` is only supported when `select` is a single block ",
           "name (e.g. \"A\"), matched against an unnamed vector or matrix ",
           "of the same shape as that block.")
    } else if (is.matrix(reference)) {
      idx <- expand.grid(row = seq_len(nrow(reference)), col = seq_len(ncol(reference)))
      ref_names <- paste0(select, "[", idx$row, ",", idx$col, "]")
    } else {
      ref_names <- paste0(select, "[", seq_along(reference), "]")
    }
    ref_df <- data.frame(parameter = ref_names, reference = as.numeric(reference))
    # left_join and restore the pre-join factor levels
    df <- dplyr::left_join(df, ref_df, by = "parameter") |>
      dplyr::mutate(parameter = factor(parameter, levels = levels(df$parameter)))
  }

  # parsing parameter levels
  df <- dplyr::mutate(df, parameter = parse_parameter(parameter))

  # random sorted subsample when there are more parameters than nshow
  pars <- unique(df$parameter)
  if (!is.null(nshow) && length(pars) > nshow) {
    pars_keep <- base::sort(sample(pars, nshow))
    df <- dplyr::filter(df, parameter %in% pars_keep)
  }

  # sort by point estimate instead of the natural parameter order
  if (sort) df <- dplyr::mutate(df, parameter = stats::reorder(parameter, m))

  # figure types
  outer_colour <- "black"
  inner_colour <- grDevices::rgb(1, 0.5, 0.1)
  if (horizontal) {
    gg <- ggplot(df, aes(y = parameter)) +
      geom_segment(aes(x = ll, xend = hh, yend = parameter), linewidth = 0.4,
                   colour = outer_colour) +
      geom_segment(aes(x = l, xend = h, yend = parameter), linewidth = 1.5,
                   colour = inner_colour) +
      geom_point(aes(x = m), size = 2, colour = outer_colour) +
      scale_y_discrete(labels = function(x) parse(text = x)) +
      labs(x = "Value", y = NULL) +
      theme_trace(legend = "none")
  } else {
    gg <- ggplot(df, aes(x = parameter)) +
      geom_segment(aes(y = ll, yend = hh, xend = parameter), linewidth = 0.4,
                   colour = outer_colour) +
      geom_segment(aes(y = l, yend = h, xend = parameter), linewidth = 1.5,
                   colour = inner_colour) +
      geom_point(aes(y = m), size = 2, colour = outer_colour) +
      scale_x_discrete(labels = function(x) parse(text = x)) +
      labs(x = NULL, y = "Value") +
      theme_trace(legend = "none")
  }

  # add reference points
  if (!is.null(reference)) {
    gg <- gg + if (horizontal) {
      geom_point(aes(x = reference), shape = 4, size = 2, colour = "forestgreen")
    } else {
      geom_point(aes(y = reference), shape = 4, size = 2, colour = "forestgreen")
    }
  }

  return(gg)
}

