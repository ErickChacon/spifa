
# CLEAN WORKSPACE AND LOAD PACKAGES --------------------------------------------

rm(list = ls())
library(datasim)
library(dplyr)
library(tidyr)
library(ggplot2)

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

data_model <- sim_model(formula = f, n = 1000, seed = 1, responses = 3)
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

output <- multi_lm(Y, X, Corr, sigmas, Cov, Beta, 1000, Sigma_proposal)

cor(output$beta)

