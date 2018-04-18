
library(tidyverse)
getwd()
Rcpp::sourceCpp("../src/correlation.cpp")
Rcpp::sourceCpp("../src/adaptive-mcmc.cpp")
source("../R/ggplot-mcmc.R")

# function_name()

k <- 25
mean <- matrix(1:k, k, 1);
N <- k * (k-1)/2
x <- rnorm(N, 0,  0.3)
chol_corr <- vec2chol_corr(x, k)
(Corr <- tcrossprod(chol_corr))
d <- sqrt(abs(rnorm(k, 2)))
D <- diag(d)
D_inv <- diag(d ^ (-1))
Sigma <- D %*% Corr %*% D


iter <- 10^6
system.time(
  samples <- adaptive_haario(mean, Sigma, iter)
)

samples_tib <- as_tibble.spmirt.list(samples, iter/2)
summary(samples_tib)
samples_long <- gather(samples_tib)

as_tibble.spmirt.list(samples, 0, 500) %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples, iter/2, 50) %>%
  gg_density_ridges(aes(fill = Parameters), scale = 6, alpha = 0.5, bandwidth = 0.25)

system.time(
  samples2 <- adaptive_haario_vanish(mean, Sigma, iter)
)

samples_tib2 <- as_tibble.spmirt.list(samples2, iter/2)
summary(samples_tib2)
samples_long2 <- gather(samples_tib2)

as_tibble.spmirt.list(samples2, 0, 500) %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples2, iter/2, 50) %>%
  gg_density_ridges(aes(fill = Parameters), scale = 6, alpha = 0.5, bandwidth = 0.25)

system.time(
  samples3 <- am_vanish_scaling(mean, Sigma, iter)
)

samples_tib3 <- as_tibble.spmirt.list(samples3, iter/2)
summary(samples_tib3)
samples_long3 <- gather(samples_tib3)

as_tibble.spmirt.list(samples3, 0, 500) %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples3, iter/2, 50) %>%
  gg_density_ridges(aes(fill = Parameters), scale = 6, alpha = 0.5, bandwidth = 0.25)


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



