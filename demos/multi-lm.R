
# CLEAN WORKSPACE AND LOAD PACKAGES --------------------------------------------

rm(list = ls())
library(datasim)
library(tidyverse)

# SIMULATE DATA ----------------------------------------------------------------

Corr <- matrix(c(1, -0.9, 0, -0.9, 1, 0.3, 0, 0.3, 1), nrow = 3)
sigmas <- c(1, 1.5, 3)
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
  geom_point(aes(col = factor(response_label)), alpha = 1 / 5)

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
source("../R/ggplot-mcmc.R")

iter <- 10^6
system.time(
  samples <- multi_lm(Y, X, iter, 0.01 * Sigma_proposal, 0.001 * Sigma_proposal)
)
samples %>% map(~ tail(.))

# apply(samples$beta, 2, mean)
# cor(samples$beta)

# Visualize traces
as_tibble(samples, 0, 100, select = "beta") %>%
  gg_trace(wrap = TRUE, alpha = 0.6)

as_tibble(samples, 0, 100, select = "beta") %>% gg_trace(alpha = 0.6)

as_tibble(samples, 0, 100, select = "corr_chol") %>% gg_trace(alpha = 0.6)

as_tibble(samples, 0, 100, select = "corr") %>% gg_trace(alpha = 0.6)

as_tibble(samples, 0, 100, select = "sigmas") %>% gg_trace(alpha = 0.6)

bla <- as_tibble(samples, iter/2, select = "sigmas")
cov(log(bla))
nrow(unique(bla)) / nrow(bla)

bla <- as_tibble(samples, iter/2, select = "corr_chol")
cov(bla)
nrow(unique(bla)) / nrow(bla)

# Visualize densities

as_tibble(samples, iter / 2, select = "corr_chol") %>%
  gg_density(aes(fill = Parameters), scale = 2, alpha = 0.5, ridges = TRUE)

as_tibble(samples, iter / 2, select = "corr") %>%
  gg_density(aes(fill = Parameters), scale = 1, alpha = 0.5, ridges = TRUE)

# Visualize credible intervals
as_tibble(samples, iter / 2, select = "beta") %>%
  summary() %>%
  mutate(param = beta) %>%
  gg_errorbarh() +
  geom_point(aes(param, Parameters), col = 3)

Corr_chol <- t(chol(Corr))
corr_chol <- Corr_chol[lower.tri(Corr_chol, diag = TRUE)]
corr <- Corr[lower.tri(Corr)]

as_tibble(samples, iter / 2, select = "corr_chol") %>%
  summary() %>%
  mutate(param = corr_chol) %>%
  gg_errorbarh() +
  geom_point(aes(param, Parameters), col = 3)

as_tibble(samples, iter / 2, select = "corr") %>%
  summary() %>%
  mutate(param = corr) %>%
  gg_errorbarh() +
  geom_point(aes(param, Parameters), col = 3)


as_tibble(samples, iter / 2 ,select = "sigmas") %>%
  summary() %>%
  mutate(param = sigmas) %>%
  gg_errorbarh() +
  geom_point(aes(param, Parameters), col = 3)


# Visualize credible intervals for all Parameters
as_tibble(samples, iter / 2) %>%
  summary() %>%
  mutate(param = c(beta, corr_chol, corr, sigmas)) %>%
  gg_errorbar() +
  geom_point(aes(Parameters, param), col = 3)

