
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
source("ggplot-mcmc.R")

iter <- 10000
samples <- multi_lm(Y, X, sigmas, iter, 0.01 * Sigma_proposal)
samples %>% map(~ tail(.))


apply(samples$beta, 2, mean)
cor(samples$beta)

# Visualize traces

samples_as_tibble(samples, "beta") %>% gg_trace(wrap = TRUE, alpha = 0.6)
samples_as_tibble(samples, "corr_chol") %>% gg_trace(alpha = 0.6)

# Visualize densities

samples_as_tibble(samples, "corr_chol") %>%
  gg_density_ridges(aes(fill = Parameters), scale = 2, alpha = 0.5)

# Visualize credible intervals
sample_summary(samples, "beta") %>%
  mutate(param = beta) %>%
  gg_errorbarh() +
  geom_point(aes(param, Parameters), col = 3)

Corr_chol <- t(chol(Corr))
corr_chol <- Corr_chol[lower.tri(Corr_chol, diag = TRUE)]

sample_summary(samples, "corr_chol") %>%
  mutate(param = corr_chol) %>%
  gg_errorbarh() +
  geom_point(aes(param, Parameters), col = 3)

# Visualize credible intervals for all Parameters
sample_summary(samples) %>%
  mutate(param = c(beta, corr_chol)) %>%
  gg_errorbar() +
  geom_point(aes(Parameters, param), col = 3)

