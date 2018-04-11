
# CLEAN WORKSPACE AND LOAD PACKAGES --------------------------------------------

rm(list = ls())
library(datasim)
library(tidyverse)

# SIMULATE DATA ----------------------------------------------------------------

Corr <- matrix(c(1, -0.9, 0, -0.9, 1, 0, 0, 0, 1), nrow = 3)
sigmas <- c(1,1,3)
D <- diag(sigmas)
Cov <- D %*% Corr %*% D
beta <- c(-10, 0, 10)

f <- list(
  mean ~ mfe(x1, beta = get("beta")) + mre(factor(id), sigma = get("Cov")),
  sd ~ I(0)
  )

data_model <- sim_model(formula = f, n = 500, responses = 3)
knitr::kable(head(data_model, 10))

# VISUALIZE DATA ---------------------------------------------------------------

ggplot(data_model, aes(x1, response)) +
  geom_point(aes(col = factor(response_label)), alpha = 1/5)

cor(matrix(data_model$mre.factor.mean, ncol = 3))

# PREPARE DATA FOR MODELLING ---------------------------------------------------

Y <- data_model %>% dplyr::select(id, response, response_label) %>%
  spread(response_label, response) %>%
  arrange(id) %>%
  dplyr::select(-id) %>%
  as.matrix()

X <- data_model %>% dplyr::select(id, matches("^x[[:digit:]]+$")) %>%
  unique() %>%
  arrange(id) %>%
  dplyr::select(-id) %>%
  as.matrix()

Beta <- matrix(beta, nrow = 1)
Sigma_proposal <- diag(1, 3)

# RUN MODEL --------------------------------------------------------------------

getwd()
Rcpp::sourceCpp("../src/multi-lm.cpp")

iter <- 5
samples <- multi_lm(Y, X, sigmas, iter, Sigma_proposal)

apply(samples$beta, 2, mean)
cor(samples$beta)

# Organize and summarise output
samples_after_burin <- samples %>% purrr::map(~ .[1:iter, ])
samples_params <- samples_after_burin
samples_params_summary <- samples_params %>%
  map(~ apply(., 2, function (x) quantile(x, c(0.025, 0.5, 0.975)))) %>%
  map(~ as_tibble(t(.))) %>%
  map(~ setNames(., make.names(names(.))))

# Visualize results
samples_params_summary$beta <- samples_params_summary$beta %>%
  mutate(param = beta)

samples_params_summary$beta %>%
  ggplot(., aes(X50., param)) +
    geom_errorbar(aes(ymin = X2.5., ymax = X97.5.)) +
    geom_point(col = 2)

