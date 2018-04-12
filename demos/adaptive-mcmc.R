
library(tidyverse)
getwd()
Rcpp::sourceCpp("../src/correlation.cpp")
Rcpp::sourceCpp("../src/adaptive-mcmc.cpp")
source("../R/ggplot-mcmc.R")

# function_name()

k <- 10
mean <- matrix(1:k, k, 1);
N <- k * (k-1)/2
x <- rnorm(N, 0,  0.3)
chol_corr <- vec2chol_corr(x, k)
(Corr <- tcrossprod(chol_corr))
d <- sqrt(abs(rnorm(k, 2)))
D <- diag(d)
D_inv <- diag(d ^ (-1))
Sigma <- D %*% Corr %*% D


iter <- 10000
system.time(
  samples <- adaptive_haario(mean, Sigma, iter)
)
class(samples) <- c("spmirt", class(samples))

as_tibble(samples) %>% gg_trace(alpha = 0.6)

# cov(samples$params) - Sigma
# apply(samples$params, 2, mean)

nrow(unique(samples$params[(iter/2+1):iter, ])) / iter * 2



