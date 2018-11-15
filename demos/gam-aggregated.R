
rm(list = ls())
setwd("../src")
getwd()

# model
library(tidyverse)
library(mgcv)

Rcpp::sourceCpp("gam.cpp")

# simulate betas
tau <- 0.01
n_knots <- 20
beta <- cumsum(cumsum(rnorm(n_knots, tau)))
plot(beta)

# simulate non-linear effect
n_ind <- 1000
data_ind <- data.frame(x = seq(-10, 10, length.out = n_ind))
sm <- smoothCon(s(x, k = n_knots, bs = "ps"), data = data_ind, knots = NULL)[[1]]
X_ind <- sm$X
D_ind <- sm$D
sigma_ind <- 1
data_ind$f <- as.numeric(X_ind %*% cbind(beta))
data_ind$y <- data_ind$f + rnorm(n_ind, sd = sigma_ind)
ggplot(data_ind, aes(x, y)) +
  geom_point() +
  geom_line(aes(x, f), col = 2, lwd = 2)

# aggregated data
n_agg <- 7
breaks <- sort(c(runif(n_agg - 1, -10, 10), -10, 10))
lengths <- diff(breaks)
data_ind$cut <- cut(data_ind$x, breaks, include.lowest = TRUE)

data_agg <- data_ind %>%
  group_by(cut) %>%
  summarise_all(mean) %>%
  extract(cut, c("low", "high"), "[(\\[](.+),(.+)[\\]]", remove = F, T)


ggplot(data_ind, aes(x, y)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(x, f), col = 4, lwd = 2) +
  geom_errorbarh(aes(xmin = low, xmax = high, y = f), data_agg, col = 4, lwd = 2) +
  geom_point(aes(x, f), data_agg, col = 4, size = 5) +
  geom_errorbarh(aes(xmin = low, xmax = high, y = y), data_agg, col = 2, lwd = 2) +
  geom_point(aes(x, y), data_agg, col = 2, size = 2)

# monte carlo approximation of the design matrix
n_mc <- 10000
data_mc <- data.frame(x = runif(n_mc, -10, 10))
data_mc$cut <- cut(data_mc$x, breaks, include.lowest = TRUE)
sm_mc <- smoothCon(s(x, k = n_knots, bs = "ps"), data = data_mc, knots = NULL)[[1]]
X_mc <- sm_mc$X
D_mc <- sm_mc$D

data_agg_mc <- data_mc %>%
  bind_cols(as.data.frame(X_mc)) %>%
  group_by(cut) %>%
  summarise_all(mean)

X_agg_mc <- dplyr::select(data_agg_mc, V1:V20) %>% as.matrix()

# run the aggregates model
bla <- gamcpp(data_agg$y, X_agg_mc, D_mc, 1, 1, niter = 1000)
beta_est <- apply(bla$beta, 2, median)
data_agg$y_est <- X_agg_mc %*% cbind(beta_est)
data_ind$y_est <- X_ind %*% cbind(beta_est)

# ggplot(data_ind, aes(x, y)) +
#   geom_point() +
#   geom_line(aes(x, f), col = 4, lwd = 1.5) +
#   geom_line(aes(x, y_est), data_ind, col = 2, size = 1.5) +
#   geom_point(aes(x, y), data_agg, col = 4, size = 3) +
#   geom_errorbarh(aes(xmin = low, xmax = high, y = y), data_agg, lwd = 2) +
#   geom_point(aes(x, y_est), data_agg, col = 2, size = 3)

bla <- ggplot(data_ind, aes(x, y)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(x, f), col = 4, lwd = 2, linetype = 2) +
  geom_errorbarh(aes(xmin = low, xmax = high, y = y), data_agg, col = 4, lwd = 1.5) +
  geom_line(aes(x, y_est), data_ind, col = 2, size = 2, linetype = 2) +
  geom_errorbarh(aes(xmin = low, xmax = high, y = y_est), data_agg, col = 2, lwd = 1.5)
bla


# pdf(file.path("gam_aggregate.pdf"), width = 7)
# print(bla)
# dev.off()




# hist(bla$sigma2 ^ 0.5, 100)
# hist(bla$tau2 ^ 0.5, 100)

# plot(bla$sigma2 ^ 0.5, type = "l")
# plot(bla$tau2 ^ 0.5, type = "l")








