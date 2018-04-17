
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


iter <- 900000
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

nrow(unique(samples_tib)) / nrow(samples_tib)
acf(samples_tib$V1, 500)

# Covariance matrix recursive update
# aux <- cbind(test[i,] - apply(test[-i, , drop = FALSE], 2, mean))
# var(test)
# variance(test)
# var(test[1:(i-1),]) * (i-2)/(i-1) + tcrossprod(aux) / i
#



