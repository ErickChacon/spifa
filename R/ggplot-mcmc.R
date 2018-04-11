
gg_trace <- function (samples, name) {
  samples %>%
    as_tibble() %>%
    setNames(paste0(name, 1:ncol(samples))) %>%
    mutate(iteration = 1:n()) %>%
    gather(varname, varvalue, -iteration) %>%
    ggplot(aes(iteration, varvalue, group = varname, col = varname)) +
      geom_path(alpha = 0.4, linetype = 1)
}


name = "bla"
samples$beta %>%
    as_tibble() %>%
    setNames(paste0(name, 1:ncol(samples$beta))) %>%
    mutate(iteration = 1:n())

params <- names(samples)
samples %>%
  purrr::map(as_tibble) %>%
  purrr::map2(params, ~ setNames(.x, paste0(.y, "[", 1:ncol(.x),"]")))

params_names <- function (name, dim, type = "vec") {
  if (type == "vec") {
    paste0(name, "[", 1:dim, "]")
  } else {
    paste0(name, "[", rep(1:dim[1], 1:dim[1]), ",", rep(1:dim[2], dim[2]:1), "]")
  }
}
params_names("beta", 8)
params_names("beta", c(3,3), "mat")
