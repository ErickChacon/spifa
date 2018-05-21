#' @title Traceplot of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' as_tibble(samples, "beta") %>% gg_trace(wrap = TRUE, alpha = 0.6)
#' as_tibble(samples, "corr_chol") %>% gg_trace(alpha = 0.6)
#'
#'
#' @export
gg_trace <- function (df, wrap = FALSE, legend = "bottom", ...) {
  df <- gather.spifa(df)
  gg <- df %>%
    ggplot(aes(iteration, Value, group = Parameters, col = Parameters)) +
      geom_path(...)
  # in splits if required
  if (wrap) {
    gg <- gg + facet_wrap(~ Parameters, ncol = 1, scales = "free")
  }
  # theme
  gg <- gg + theme(legend.position = legend)
  return(gg)
}

#' @title Densities of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' as_tibble(samples, "corr_chol") %>%
#'   gg_density_ridges(aes(fill = Parameters), scale = 2, alpha = 0.5)
#'
#' @export
gg_density <- function (df, ..., ridges = FALSE) {
  df <- gather.spifa(df)
  df <- df %>%
    group_by(Parameters) %>%
    mutate(median = quantile(Value, 0.5))
  if (ridges) {
    gg <- df %>%
      ggplot(aes(Value, Parameters, group = Parameters)) +
        ggridges::geom_density_ridges(...)
  } else  {
    gg <- df %>%
      ggplot(aes(Value, fill = Parameters)) +
      geom_density(...) +
      facet_wrap(~ Parameters, scales = "free")
  }
  # theme
  gg <- gg + theme(legend.position = "none")
  return(gg)
}

#' @title 2D Densities of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' as_tibble(samples, "corr_chol") %>%
#'   gg_density_ridges(aes(fill = Parameters), scale = 2, alpha = 0.5)
#'
#' @export
gg_density2d <- function (samples, var1, var2, each = NULL,
                          keys = c("group", "Parameter"), highlight = NULL,
                          ncol = NULL, ...) {

  if (!is.null(highlight)) {
    aux_samples <- samples[1,]
    aux_samples[1, ] <- highlight
  }

  if (!is.null(each)) {
    samples <- gather.spifa(samples, each, keys)
    aux_samples <- gather.spifa(aux_samples, each, keys)
  }

  gg <- ggplot(samples, aes_(substitute(var1), substitute(var2))) +
    stat_density2d(aes(fill = log(..level..)),
                   geom = 'polygon', col = "black", ...) +
    # scale_fill_continuous(low="green",high="red") +
    guides(alpha="none")
  # +
  #   geom_point(...)

  if (!is.null(highlight)) {
    gg <- gg + geom_point(data = aux_samples, col = 2)
  }

  if (!is.null(each)) {
    if (!is.null(ncol)) {
      gg <- gg + facet_wrap(~ groups, scales = "free", ncol = ncol)
    } else {
      gg <- gg + facet_wrap(~ groups, scales = "free")
    }
  }


  return(gg)
}



#' @title 2D Scatterplot of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#'
#' @export
gg_scatter <- function (samples, var1, var2, each = NULL,
                        keys = c("group", "Parameter"), highlight = NULL,
                        ncol = NULL, points_alpha = 0.5) {

  if (!is.null(highlight)) {
    aux_samples <- samples[1,]
    aux_samples[1, ] <- highlight
  }

  if (!is.null(each)) {
    samples <- gather.spifa(samples, each, keys)
    aux_samples <- gather.spifa(aux_samples, each, keys)
  }

  gg <- ggplot(samples, aes_(substitute(var1), substitute(var2))) +
    geom_point(alpha = points_alpha) +
    geom_path(alpha = 0.4, linetype = 2)

  if (!is.null(highlight)) {
    gg <- gg + geom_point(data = aux_samples, col = 2, size = 2)
  }

  if (!is.null(each)) {
    if (!is.null(ncol)) {
      gg <- gg + facet_wrap(~ groups, scales = "free", ncol = ncol)
    } else {
      gg <- gg + facet_wrap(~ groups, scales = "free")
    }
  }

  return(gg)
}






#' @title Horizontal Errorbar Plot of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' summary(samples, "corr_chol") %>%
#'   mutate(param = corr_chol) %>%
#'   gg_errorbarh() +
#'   geom_point(aes(param, Parameters), col = 3)
#'
#' @export
gg_errorbarh <- function (df_summary, sorted = FALSE,
                          colors = c(rgb(1,0.5,0.1), "black"), ...) {

  if (sorted) {
    gg <- df_summary %>% ggplot(., aes(`50%`, `50%`))
  } else {
    gg <- df_summary %>% ggplot(., aes(`50%`, Parameters))
  }

  gg <- gg +
    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, col = "95%"),
                   height = 0, ...) +
    geom_errorbarh(aes(xmin = `10%`, xmax = `90%`, col = "80%"), size = 2,
                   height = 0, ...) +
    geom_point(size = 2)
  # colors
  gg <- gg + scale_colour_manual(values = colors)
  # labels
  gg <- gg + labs(colour = "Credible Intervals:", x = "Value")
  # theme
  gg <- gg + theme(legend.position = "bottom")
  return(gg)
}

#' @title Errorbar Plot of Samples
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' summary(samples, "corr_chol") %>%
#'   mutate(param = corr_chol) %>%
#'   gg_errorbar() +
#'   geom_point(aes(Parameters, param), col = 3)
#'
#' @export
gg_errorbar <- function (df_summary, sorted = TRUE,
                         colors = c(rgb(1,0.5,0.1), "black"), ...) {

  if (sorted) {
    gg <- df_summary %>% ggplot(., aes(`50%`, `50%`))
  } else {
    gg <- df_summary %>% ggplot(., aes(Parameters, `50%`))
  }

  gg <- gg +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, col = "95%"),
                   width = 0, ...) +
    geom_errorbar(aes(ymin = `10%`, ymax = `90%`, col = "80%"), size = 2,
                   width = 0, ...) +
    geom_point(size = 2)
  # colors
  gg <- gg + scale_colour_manual(values = colors)
  # labels
  gg <- gg + labs(colour = "Credible Intervals:", x = "Value")
  # theme
  gg <- gg + theme(legend.position = "bottom")
  return(gg)
}
