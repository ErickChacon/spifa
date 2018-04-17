
library(tidyverse)
getwd()
Rcpp::sourceCpp("../src/correlation.cpp")
Rcpp::sourceCpp("../src/adaptive-mcmc.cpp")
source("../R/ggplot-mcmc.R")

# function_name()

k <- 15
mean <- matrix(1:k, k, 1);
N <- k * (k-1)/2
x <- rnorm(N, 0,  0.3)
chol_corr <- vec2chol_corr(x, k)
(Corr <- tcrossprod(chol_corr))
d <- sqrt(abs(rnorm(k, 2)))
D <- diag(d)
D_inv <- diag(d ^ (-1))
Sigma <- D %*% Corr %*% D


iter <- 500000
system.time(
  samples <- adaptive_haario(mean, Sigma, iter)
)
class(samples) <- c("spmirt", class(samples))

samples_tib <- as_tibble(samples)
# samples_tib %>% gg_trace(alpha = 0.6)

samples_tib[(iter/3+1):iter, ] %>%
  gg_density_ridges(aes(fill = Parameters), scale = 6, alpha = 0.5, bandwidth = 0.25)

# hist(samples_tib$V1[(iter/10+1):iter], 100, probability = TRUE)
# lines(density(samples_tib$V1[(iter/10+1):iter], bw = 0.25))
# hist(rnorm(iter, 1, 1), 100)
# cov(samples$params) - Sigma
# apply(samples$params, 2, mean)

# nrow(samples)

nrow(unique(samples$params[(iter/2+1):iter, ])) / iter * 2
nrow(unique(samples$params)) / iter

# n = 8
# i = n/2
# test <- matrix(rnorm(n), ncol = 2)
#
# aux <- cbind(test[i,] - apply(test[-i, , drop = FALSE], 2, mean))
# var(test)
# variance(test)
# var(test[1:(i-1),]) * (i-2)/(i-1) + tcrossprod(aux) / i
#



