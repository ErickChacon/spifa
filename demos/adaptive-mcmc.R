
library(tidyverse)
getwd()
Rcpp::sourceCpp("../src/correlation.cpp")
Rcpp::sourceCpp("../src/adaptive-mcmc.cpp")
source("../R/ggplot-mcmc.R")

k <- 35
mean <- matrix(1:k, k, 1);
N <- k * (k - 1) / 2
x <- rnorm(N, 0,  0.3)
chol_corr <- vec2chol_corr(x, k)
(Corr <- tcrossprod(chol_corr))
d <- sqrt(abs(rnorm(k, 2)))
D <- diag(d)
D_inv <- diag(d ^ (-1))
Sigma <- D %*% Corr %*% D

iter <- 1.5 * 10 ^ 6
system.time(
  samples <- adaptive_haario(mean, Sigma, iter)
)

samples_tib <- as_tibble.spmirt.list(samples, iter / 2)
summary(samples_tib)
samples_long <- gather(samples_tib)

as_tibble.spmirt.list(samples, 0, 500) %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples, iter / 2, 50) %>%
  gg_density(aes(fill = Parameters), scale = 6, alpha = 0.5, bandwidth = 0.25,
             ridges = TRUE)

system.time(
  samples2 <- adaptive_haario_vanish(mean, Sigma, iter)
)

samples_tib2 <- as_tibble.spmirt.list(samples2, iter/2)
summary(samples_tib2)
samples_long2 <- gather(samples_tib2)

as_tibble.spmirt.list(samples2, 0, 500) %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples2, iter/2, 50) %>%
  gg_density(aes(fill = Parameters), scale = 6, alpha = 0.5, bandwidth = 0.25,
             ridges = TRUE)

system.time(
  samples3 <- am_vanish_scaling(mean, Sigma, iter)
)

samples_tib3 <- as_tibble.spmirt.list(samples3, iter/2)
summary(samples_tib3)
samples_long3 <- gather(samples_tib3)

as_tibble.spmirt.list(samples3, 0, 500) %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples3, iter/2, 50) %>%
  gg_density(aes(fill = Parameters), scale = 6, alpha = 0.5, bandwidth = 0.25,
             ridges = TRUE)


nrow(unique(samples_tib)) / nrow(samples_tib)
nrow(unique(samples_tib2)) / nrow(samples_tib2)
nrow(unique(samples_tib3)) / nrow(samples_tib3)

coda::effectiveSize(samples_tib)
coda::effectiveSize(samples_tib2)
coda::effectiveSize(samples_tib3)

# cov(samples_tib) - Sigma
# cov(samples_tib2) - Sigma

# nrow(unique(samples_tib)) / nrow(samples_tib)
# acf(samples_tib$V3, 500)
# acf(samples_tib2$V3, 500)
# acf(samples_tib3$V3, 500)
# acf(samples_tib$V2, 500)
# acf(samples_tib2$V2, 500)

# Covariance matrix recursive update
# aux <- cbind(test[i,] - apply(test[-i, , drop = FALSE], 2, mean))
# var(test)
# variance(test)
# var(test[1:(i-1),]) * (i-2)/(i-1) + tcrossprod(aux) / i
#

# Analizing gamma series for adaptive mcmc
# i <- 0:1000
# C <- 1
# alpha <- 1
# plot(i, C / (i+1)^alpha, type = "l")
# C <- 0.9
# alpha <- 1
# lines(i, C / (i+1)^alpha, col = 2)
# C <- 0.9
# alpha <- 0.9
# lines(i, C / (i+1)^alpha, col = 3)
# C <- 0.8
# alpha <- 0.5
# lines(i, C / (i+1)^alpha, col = 4)
#
#
#
# C <- 1
# alpha <- 1
# i <- 1:10000
# plot(i, purrr::map_dbl(i, ~ sum(C / (1:.)^ alpha)), type = "l", col = 2)
# lines(i, purrr::map_dbl(i, ~ sum((C / (1:.)^ alpha) ^ (1+1))))
#
# C <- 0.9
# alpha <- 1
# # i <- 1:10000
# # plot(i, purrr::map_dbl(i, ~ sum(C / (1:.)^ alpha)), type = "l", col = 2)
# # lines(i, purrr::map_dbl(i, ~ sum((C / (1:.)^ alpha) ^ (1+1))))
# #
# # C <- 0.9
# # alpha <- 0.9
# # i <- 1:10000
# # plot(i, purrr::map_dbl(i, ~ sum(C / (1:.)^ alpha)), type = "l", col = 2)
# # lines(i, purrr::map_dbl(i, ~ sum((C / (1:.)^ alpha) ^ (1+1))))
# #
# # C <- 0.9
# # alpha <- 0.5
# # i <- 1:10000
# # plot(i, purrr::map_dbl(i, ~ sum(C / (1:.)^ alpha)), type = "l", col = 2)
# # lines(i, purrr::map_dbl(i, ~ sum((C / (1:.)^ alpha) ^ (1+1))))
