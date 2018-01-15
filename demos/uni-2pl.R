# CLEAN WORKSPACE AND LOAD PACKAGES --------------------------------------------

rm(list = ls())

library(tidyverse)
library(day2day)
# library(spmirt)
library(coda)
# library(rstan)
# library(Matrix)

# SIMULATE DATA ----------------------------------------------------------------

q <- 10
n <- 300
nq <- n * q
sigma <- 1
discrimination <- seq(0.4, 1.5, length.out = q)
# discrimination <- rep(1, q)
difficulty <- (1:q - 5)/10 * 2
# difficulty <- rep(0, q)
ability <- rnorm(n, 0, sigma)
data_long <- expand.grid(subject = 1:n, item = 1:q) %>% mutate(
  ability = ability[subject],
  difficulty = difficulty[item],
  discr = discrimination[item],
  y = rbinom(n * q, 1, psych::logistic(1.702 * (discr * ability - difficulty))),
  # y = rbinom(n * q, 1, psych::logistic(discr * (ability - difficulty))),
  # y_logit = rbinom(n * q, 1, psych::logistic(discr * ability - difficulty)),
  # y_probit = rbinom(n * q, 1, psych::logistic(1.702 * (discr * ability - difficulty)))
  )

bla <- data_long %>% as_tibble() %>% group_by(subject) %>%
  summarize(a1 = mean(y), abil = mean(ability))
plot(bla$abil, bla$a1)

data_long %>% as_tibble() %>% group_by(item) %>% summarize(easiness = mean(y))
bla <- data_long %>% as_tibble() %>% group_by(subject) %>% summarize(ability = mean(y))
table(bla$ability)

# FIT MODEL WITH SPMIRT --------------------------------------------------------


ini <- Sys.time()
samples <- ifa_gibbs(data_long$y, n, q, 5000)
print(Sys.time() - ini)
# 45

# ini <- Sys.time()
# samples <- ifa_gibbs_no(data_long$y, n, q, 5000)
# print(Sys.time() - ini)
#

# samples2 <- samples %>% purrr::map(~ .[1001:2000, ])
samples2 <- samples %>% purrr::map(~ .[2501:5000, ])
means <- purrr::map(samples2, ~ apply(., 2, mean))
medians <- purrr::map(samples2, ~ apply(., 2, median))

plot(difficulty, -means$c)
points(difficulty, -medians$c, col = 2)
abline(0, 1)

c_sam <- mcmc(-samples2$c)
plot(c_sam)

plot(discrimination, means$a)
abline(0, 1)
abline(0, -1)

a_sam <- mcmc(samples2$a)
plot(a_sam)

plot(ability, means$theta)
abline(0, 1)
abline(0, -1)

# theta_sam <- mcmc(samples2$theta)
# plot(theta_sam)
# plot(samples$c[, 1], type = "l")
# hist(samples$c[, 1])

plot(samples2$z[,1], samples2$theta[,1])

plot(means$z, data_long$y)
abil <- tapply(means$z, data_long$subject, mean)
plot(ability, abil)
abline(0, 1)


# psych::logistic(1.702)
# pnorm(1.702)


# x <- seq(-3, 3, 0.1)
# plot(psych::logistic(1.702 * x), pnorm(x))
# abline(0, 1)
