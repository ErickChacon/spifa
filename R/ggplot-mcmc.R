
# as_tibble <- function (x, ...) UseMethod("as_tibble", x)

#' @title Convert Samples to Tibble
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
#' as_tibble(samples)
#' as_tibble(samples, "beta")
#'
#' @export
as_tibble.spmirt <- function (samples, select = names(samples)) {
  samples <- samples[select]
  params <- names(samples)
  samples %>%
    purrr::map(as_tibble) %>%
    dplyr::bind_cols()
}

#' @title Summary of Samples
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
#' summary(samples)
#'
#' @export
summary.spmirt <- function (samples, select = names(samples)) {

  df <- samples_as_tibble(samples, select)
  df_names <- names(df)
  df <- df %>%
    map(function (x) quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))) %>%
    reduce(rbind) %>%
    as_tibble() %>%
    mutate(Parameters = df_names)
  df <- df[c(ncol(df), 1:(ncol(df)-1))]
  return(df)
}

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
#' samples_as_tibble(samples, "beta") %>% gg_trace(wrap = TRUE, alpha = 0.6)
#' samples_as_tibble(samples, "corr_chol") %>% gg_trace(alpha = 0.6)
#'
#'
#' @export
gg_trace <- function (df, wrap = FALSE, ...) {
  gg <- df %>%
    as_tibble() %>%
    mutate(iteration = 1:n()) %>%
    gather(Parameters, Value, -iteration) %>%
    ggplot(aes(iteration, Value, group = Parameters, col = Parameters)) +
      geom_path(...)
  # in splits if required
  if (wrap) {
    gg <- gg + facet_wrap(~ Parameters, ncol = 1, scales = "free")
  }
  # theme
  gg <- gg + theme(legend.position = "bottom")
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
#' samples_as_tibble(samples, "corr_chol") %>%
#'   gg_density_ridges(aes(fill = Parameters), scale = 2, alpha = 0.5)
#'
#' @export
gg_density_ridges <- function (df, ...) {
  gg <- df %>%
    as_tibble() %>%
    mutate(iteration = 1:n()) %>%
    gather(Parameters, Value, -iteration) %>%
    ggplot(aes(Value, Parameters, group = Parameters)) +
      ggridges::geom_density_ridges(...)
  # theme
  gg <- gg + theme(legend.position = "none")
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
#'   mutate(param = blaa) %>%
#'   gg_errorbarh() +
#'   geom_point(aes(param, Parameters), col = 3)
#'
#' @export
gg_errorbarh <- function (df_summary, colors = c(rgb(1,0.5,0.1), "black"), ...) {
  gg <- df_summary %>%
    ggplot(., aes(`50%`, Parameters)) +
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
gg_errorbar <- function (df_summary, colors = c(rgb(1,0.5,0.1), "black"), ...) {
  gg <- df_summary %>%
    ggplot(., aes(Parameters, `50%`)) +
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
