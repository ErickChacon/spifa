# CLEAN WORKSPACE AND LOAD PACKAGES --------------------------------------------

rm(list = ls())

library(spmirt)
library(ggplot2)
library(datasim)
library(tidyverse)

# SIMULATE DATA ----------------------------------------------------------------

n <- 300

f <- list(
  prob ~ I(0) +
    gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.04),
       sigma2 = 1),
  size ~ I(1)
  )

data <- sim_model(formula = f,
                  link_inv = list(pnorm, identity),
                  generator = rbinom,
                  n = n,
                  seed = 2
                  )
data <- dplyr::rename(data, gp = gp.list.prob)

ggplot(data, aes(s1, s2)) +
  geom_point(aes(col = gp, size = gp)) +
  scale_colour_distiller(palette = "RdYlBu")

ggplot(data, aes(s1, s2)) +
  geom_point(aes(col = factor(response)), size = 2)

data$gp2 = data$gp + rnorm(300)


vg <- gstat::variogram(gp2 ~ 1, ~ s1 + s2, data, cutoff = 0.7, width = 0.005)
ggplot(vg, aes(dist, gamma, weight = np)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.7))

vg <- gstat::variogram(response ~ 1, ~ s1 + s2, data, cutoff = 0.7, width = 0.005)
ggplot(vg, aes(dist, gamma)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.7))

noise <- rnorm(n)

##        V1        V2
##  433.6714 1068.1642
##        V1        V2
##  423.8379 1083.1148
##        V1        V2
##  684.2511 1249.9670
# RUN MODELS -------------------------------------------------------------------

Rcpp::sourceCpp("../src/gp-gibss-adap.cpp")
source("../R/ggplot-mcmc.R")

# 50000 -> 17 minutes

# First naive adaptive
iter <- 5 * 10 ^ 4
dist <- as.matrix(dist(dplyr::select(data, s1, s2)))
sigma_prop <- matrix(c(0.138, -0.023, -0.023, 0.1), 2) * 2.38 ^ 2 / 2

system.time(
  samples <- probit_gp_adap(data$response, dist, c(log(1), log(0.04)), iter,
                            sigma_prop)
)

samples_tib <- as_tibble.spmirt.list(samples, iter/2, select = "param")
summary(samples_tib)
samples_long <- gather(samples_tib)

as_tibble.spmirt.list(samples, 0, 1, "param") %>%
  gg_trace(alpha = 0.6, wrap = TRUE)

as_tibble.spmirt.list(samples, iter/2, 1, "param") %>%
  gg_density(alpha = 0.5)

as_tibble.spmirt.list(samples, iter/2, 15, "param") %>%
  gg_density2d(V1, V2)

as_tibble.spmirt.list(samples, iter/2, 15, "param") %>%
  gg_density2d(log(V1), log(V2))

nrow(unique(samples_tib)) / nrow(samples_tib)

coda::effectiveSize(log(samples_tib))
acf(samples_tib)

cov(log(samples_tib))

# Adaptive with stochastic approximation series

system.time(
  samples2 <- probit_gp_am(data$response, dist, c(log(1), log(0.04)), iter,
                            sigma_prop)
)

samples_tib2 <- as_tibble.spmirt.list(samples2, iter/2, select = "param")
summary(samples_tib2)
samples_long2 <- gather(samples_tib2)

as_tibble.spmirt.list(samples2, 0, 1, "param") %>%
  gg_trace(alpha = 0.6, wrap = TRUE)

as_tibble.spmirt.list(samples2, iter/2, 1, "param") %>%
  gg_density(alpha = 0.5)

as_tibble.spmirt.list(samples2, iter/2, 15, "param") %>%
  gg_density2d(V1, V2)

as_tibble.spmirt.list(samples2, iter/2, 15, "param") %>%
  gg_density2d(log(V1), log(V2))

nrow(unique(samples_tib2)) / nrow(samples_tib2)

coda::effectiveSize(log(samples_tib2))
acf(samples_tib2)

cov(log(samples_tib2))

# Adaptive with stochastic approximation series

system.time(
  samples3 <- probit_gp_am_scale(data$response, dist, c(log(1), log(0.04)), iter,
                            sigma_prop, 0.35)
)

samples_tib3 <- as_tibble.spmirt.list(samples3, iter/2, select = "param")
summary(samples_tib3)
samples_long3 <- gather(samples_tib3)

as_tibble.spmirt.list(samples3, 0, 1, "param") %>%
  gg_trace(alpha = 0.6, wrap = TRUE)

as_tibble.spmirt.list(samples3, iter/2, 1, "param") %>%
  gg_density(alpha = 0.5)

as_tibble.spmirt.list(samples3, iter/2, 15, "param") %>%
  gg_density2d(V1, V2)

as_tibble.spmirt.list(samples3, iter/2, 15, "param") %>%
  gg_density2d(log(V1), log(V2))

nrow(unique(samples_tib3)) / nrow(samples_tib3)

coda::effectiveSize(log(samples_tib3))
acf(samples_tib3)

cov(log(samples_tib3))
#        V1        V2
#  671.1306 1195.5363




# system.time(
#   out <- probit_gp_chol2(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
# )
# system.time(
#   out0 <- probit_gp(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
# )

# plot(out$param[, 2])
# abline(h = 0.03, col = 2)



# # Compare priors with posteriors
# # phi
# hist(exp(rnorm(50000, log(0.03), 0.4)), 200, freq = FALSE )
#

# system.time(
#   out <- probit_gp_chol2(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
# )
# system.time(
#   out0 <- probit_gp(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
# )

# plot(out$param[, 2])
# abline(h = 0.03, col = 2)



# # Compare priors with posteriors
# # phi
# hist(exp(rnorm(50000, log(0.03), 0.4)), 200, freq = FALSE )
# lines(density(out$param[(iter/2):iter, 2]), col = 2)
# quantile(out$param[(iter/2):iter, 2], c(0.025, 0.975))
#
# # sigma^2
# hist(exp(rnorm(50000, log(1), 0.4)), 200, freq = FALSE, xlim = c(0, 5))
# # hist(exp(rnorm(5000, log(1), 0.3)), 100, freq = FALSE)
# lines(density(out$param[(iter/2):iter, 1]), col = 2)
# quantile(out$param[(iter/2):iter, 1], c(0.025, 0.975))
#
# # phi
# hist(exp(rnorm(50000, log(0.03), 0.4)), 200, freq = FALSE )
# lines(density(out$param[(iter/2):iter, 2]), col = 2)
# quantile(out$param[(iter/2):iter, 2], c(0.025, 0.975))
#
# # sigma^2
# hist(exp(rnorm(50000, log(1), 0.4)), 200, freq = FALSE, xlim = c(0, 5))
# # hist(exp(rnorm(5000, log(1), 0.3)), 100, freq = FALSE)
# lines(density(out$param[(iter/2):iter, 1]), col = 2)
# quantile(out$param[(iter/2):iter, 1], c(0.025, 0.975))
#
# acf(out$param[(iter/2):iter, 1])
# acf(out$param[(iter/2):iter, 2])
#
# acf(out$param[seq((iter/2), iter, 10), 1])
# # acf(out$param[seq((iter/2), iter, 10), 2])
#
